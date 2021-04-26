

%%%% EDIT: NEW APPROACH : SHOW AIRS AND MODEL AS AIRS ON THE SAME AIRS
%%%% GRID USING THE MOSAIC APPROACH TO SHOW THE AIRS FOOTPRINTS

%%%% This is called the "airs_grid" of "AG" version

% making a new model as airs processor...

% %%
% modelday = datenum(2013,07,15);
% yyyymmdd = '20150705';
% airsdirec = '/Volumes/SDBlue/data/sgwex/AIRS/three_corner_overlapping_granule_pairs/matlab/';
% modeldirec = '/Volumes/BLUE/data/sgwex/Model/JustT/';
% savedirec = '~/Desktop/';


function FUN_sgwex_model_as_airs_regrid_and_3dst_AG(yyyymmdd,airsdirec,modeldirec,savedirec)


% for each hour today, calculate the nearest airs overpass:
airsfiles = dir([airsdirec '*_AirsSG.mat']);
airsgrans = cat(1,airsfiles(:).name);
airsgrans = cellstr(airsgrans(:,1:13));
airstimes = datenum(airsgrans,'yyyymmdd-HHMM');
loadedairsgranule = '-----';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load and AIRS granule to get a retrieval noise profile from.
% 20150620-1700 looks good.
% best to use a real granule cos the noise is different for each altitude.
A = load([airsdirec '20150620-1700_AirsSG.mat']);
[A.Airs.tp,~] = nph_airs_4dp_detrend(A.Airs.T,1,4);
airsnoise = A.Airs.tp(:,136:270,:);
airsnoise = cat(2,airsnoise,airsnoise);
% randomise each layer to get rid of the join:
asz = size(A.Airs.T);
for z = 1:asz(3)
    layer = airsnoise(:,:,z);
    layer = reshape(layer(randperm(prod(asz(1:2)))),asz(1:2));
    airsnoise(:,:,z) = layer;
end
% and you're done!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Loading Model Lat Lon Alt grid...')
% saves time regridding and loading...
M = load([modeldirec 'ModelLatLonGrid.mat']);
ModelLatLonGrid = M.ModelGrid;
clear M

%%%% FOR EACH HOUR, LOAD HOURLY MODEL TIMESTEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modelday = datenum(yyyymmdd,'yyyymmdd');

for h = 0:23
    
    granule = [yyyymmdd '-' sprintf('%02d',h) '00'];
    
    disp(' ')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(granule)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% 1 - load full res model
    %     load([modeldirec granule '_SG_uvwtpd.mat']);
    try
%         disp('- Loading Model...')
%         M = load([modeldirec granule '_Model_3DST.mat']);
%         Model = M.Model;
%         clear M
        disp('- Loading Model Temperature field...')
        M = load([modeldirec granule '_Model_T.mat']);
        Model = M.Model;
        Model.Lat = ModelLatLonGrid.Lat;
        Model.Lon = ModelLatLonGrid.Lon;
        Model.Alt = ModelLatLonGrid.Alt;
        clear M
    catch err
        disp(['**** Error loading Model for ' granule ' ****'])
        disp(err.identifier)
        continue
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% 0 - before anything else, find the nearest airs time,
    %%%% and load AIRS granule:
    
    [~,ind] = min(abs(datenum(Model.Time) - airstimes));
    nearestgranule = airsfiles(ind).name(1:13);
    
    % check if this is the one you already have loaded:
    if ~strcmpi(loadedairsgranule,nearestgranule)
        disp(['- Loading NEW nearest AIRS granule ' nearestgranule])
        A = load([airsdirec nearestgranule '_AirsSG.mat']);
        Airs = A.Airs;
        clear A
        loadedairsgranule = nearestgranule;
    else
        % what airs is it?
        disp(['- Using nearest AIRS granule ' nearestgranule ])
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% 2. simulate the airs horz footrpint in the model
    
    horz_fwhm = [13.5 13.5]; % km
    horz_grid_spacing = [1.5 1.5]; % km
    sig = (horz_fwhm ./ horz_grid_spacing) ./ 2.355;
    Model.T_sm = double(imgaussfilt(Model.T,sig,...
        'padding','replicate',...
        'filterdomain','spatial'));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. regrid model temp to the AIRS scantrack.
%     disp('- Regridding model to airs scan track...')
    
    clon = Model.Lon(300,400);
    clat = Model.Lat(300,400);
    
    % AIRS distance from model centre
    [Airs.d,Airs.az] = distance(clat,clon,Airs.Lat,Airs.Lon);
    Airs.xdist = deg2km(Airs.d).*sind(Airs.az);
    Airs.ydist = deg2km(Airs.d).*cosd(Airs.az);
    
    % Now you only need to make a GRIDDED interpolant for the model:
    mx = linspace(-600,598.5,800);
    my = linspace(-450,448.5,600);
    mz = Model.Alt;
    Fh_model = griddedInterpolant({my,mx,mz},Model.T_sm,'linear','nearest');
    % ^ note the nearest neighbour extrapolation! this makes sense because
    % the model edges are basically the low res global forecast.
    
    % finally, regrid to the AIRS scantrack
    szm = size(Model.T);
    sza = size(Airs.T);
    X = repmat(Airs.xdist,1,1,szm(3));
    Y = repmat(Airs.ydist,1,1,szm(3));
    
    % evaluatue on the airs scan track:
    Model.Ta = nan(sza(1),sza(2),szm(3));
    Z = permute(repmat(Model.Alt(:),1,sza(1),sza(2)),[2 3 1]);
    Model.Ta = Fh_model(Y,X,Z);
    
% % % % % % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % %     %%%% 3.75 WEIGHT PIXEL-SCALE VARIATIONS TO NEAR THE ISLAND
% % % % % % % % %     % in consistency with my latest AIRS analysis, apply a weighted average
% % % % % % % % %     % between the raw data and a smoothed one with increasing distance from
% % % % % % % % %     % the island.
% % % % % % % % %     
% % % % % % % % %     % make a smoothed version:
% % % % % % % % %     smoo = smoothn(Model.Ta,[3 3 1]);
% % % % % % % % %     
% % % % % % % % %     % Average the smoothed and non-smoothed temperatures, but use range
% % % % % % % % %     % from the island as a weighting for whether we want more of the
% % % % % % % % %     % smoothed or more of the unsmoothed.
% % % % % % % % %     rng = quadadd(X-100,Y);
% % % % % % % % %     rng(rng > 400) = 400; % set a limit for transitioning to the smoothed.
% % % % % % % % %     rng = rng.^0.5; % get a bit faster roll off
% % % % % % % % %     rng = rng ./ max(rng(:));
% % % % % % % % %     
% % % % % % % % %     Model.Ta = (Model.Ta.*(1-rng)) + (smoo.*(rng));
% % % % % % % % %     % each location is esentially a two-element weighted average between
% % % % % % % % %     % the smoothed T and the raw T at that location.
% % % % % % % % %     % If the weights add up to 1, then you just multiply each value by the
% % % % % % % % %     % weighting and add it up.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4. Apply the airs vertical resolution. This should be faster now because
    % 90x270x50 is a lot smaller than 600x800x118!
%     disp('- Applying airs vertical resolution(s)...')
    
    vert_resi = floor(airs3d_vert_res(Model.Alt,'nearest')); % give benefit of the doubt ;)
    vert_grid_spacing = nanmean(diff(Model.Alt));
    uvr = unique(vert_resi);
    Model.Ta_sm = zeros(size(Model.Ta));
    
    % substitute any nans resulting from the airs regridding:
    nanlocs = isnan(Model.Ta);
    for z = 1:szm(3)
        layer = Model.Ta(:,:,z);
        layer(isnan(layer)) = nanmedian(linearise(layer));
        Model.Ta(:,:,z) = layer;
    end
    
    % now for each vertical resolution:
    for i = 1:length(uvr)
        
%         disp(uvr(i))
        
        % compute 3d fwhm and grid spacing:
        fwhm = [0.1 0.1 uvr(i)];
        grid_spacing = [1 1 vert_grid_spacing];
        
        % smooth the model for this vert_res:
        sig = (fwhm ./ grid_spacing) ./ 2.355;
        fsize = 2*ceil(2*sig)+1;
        fsize(1:2) = 1;
        smoo = imgaussfilt3(Model.Ta,sig,...
            'padding','symmetric',...
            'filterdomain','spatial');
        
        % identify all the altitudes associated with this vert res:
        zinds = find(vert_resi == uvr(i));
        
        % and assign the smoothed values:
        Model.Ta_sm(:,:,zinds) = smoo(:,:,zinds);
        
    end
    
    % put back the nans for consistency:
    Model.Ta_sm(nanlocs) = NaN;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% 4.5 SIMULATE AIRS RETRIEVAL NOISE
    
    %%%% note that this must be done AFTER the vertical resolution has been
    %%%% applied for compatability with AIRS.
    % before, I was doing it before the vert res, so all the noise what
    % getting smoothed out.
    %%%%%%%%%%%
    % ok so if we're gonna do area averaged 3DST values, we need to simulate
    % some AIRS retrieval noise in the model-as-AIRS. If you can't beat
    % them, join them!
    
    % Also use the 3x3 boxcar to smooth. Let's not go crazy.
    airsnoise_sm = smoothn(airsnoise,[3 3 1]);
    
    % put the airs noise on the same grid as the model:
    Fnoise = griddedInterpolant({1:sza(1),1:sza(2),Airs.Alt},airsnoise_sm,'linear','nearest');
    
    % and add it on
    Model.Ta_sm = Model.Ta_sm + Fnoise({1:sza(1),1:sza(2),Model.Alt});
    
    % now, for consistency with the AIRS approach, make a smoothed version
    % and average it with the raw version using a weighting based on
    % distance from the island.
    smoo = smoothn(Model.Ta_sm,[3 3 1]);
    rng = repmat(quadadd(Airs.xdist-100,Airs.ydist),1,1,50);
    rng(rng > 400) = 400; % set a limit for transitioning to the smoothed.
    rng = rng.^0.5; % get a bit faster roll off
    rng = rng ./ max(rng(:));
    
    Model.Ta_sm = (Model.Ta_sm.*(1-rng)) + (smoo.*(rng));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% 5. Extract GW perturbations IN CROSS TRACK DIRECTION
    %     disp('- Extracting GW perturbations...')
    
    % note the change in order. by using nearest neighbour extrapolation in the
    % regridding to the airs scan track, we can use the full cross-track
    % detrending.
    
    % i was finding that when the AIRS is detrended on the regular grid, we get
    % a lot of weird artefacts at the grid edge. by using the full XT width,
    % these should be reduce as we're only taking the central region.
    
    [Model.Tpg,Model.bga] = nph_airs_4dp_detrend(Model.Ta_sm,1,4);
    % ^ only need the background, plus this will regrid easier in the next
    % step.
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% 6. Prepare a regular grid for 3DST
    % this needs changing a bit from previous versions to allow for shorter
    % LH waves.
%     disp('- Putting model as airs on regular grid...')
    
    % define a new grid size
    gsz = [90 120 length(Model.Alt)];
    g1 = linspace(-450,450,gsz(1));
    g2 = linspace(-600,600,gsz(2));
    
    % and find the lat/lon of it:
    [G1,G2] = ndgrid(g1,g2);
    d = km2deg(quadadd(G1,G2));
    az = atan2d(G2,G1);
    [Model.Latg,Model.Long] = reckon(clat,clon,d,az);
    
    % get airs scan track interpolant: (have to use scattered interpolant
    % i'm afraid...)
    Fhorz = scatteredInterpolant(Airs.ydist(:),Airs.xdist(:),zeros(size(Airs.xdist(:))),'linear','none');
    
    Model.Tg = nan(gsz);
    Model.bgg = nan(gsz);
    for z = 1:length(Model.Alt)
        % raw temperature
        Fhorz.Values = linearise(Model.Ta_sm(:,:,z));
        Model.Tg(:,:,z) = Fhorz({g1,g2});
        % fitted background
        Fhorz.Values = linearise(Model.bga(:,:,z));
        Model.bgg(:,:,z) = Fhorz({g1,g2});
    end
    
    % finally, subtract to get gridded temperature perturbations:
    Model.Tpg = Model.Tg - Model.bgg;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% 7. For compatibility with airs, remove altitudes where we don't
    %%%% trust the airs retrieval.
    
%     disp('- Trimming to useable height range...')
    
    lower_height_index = 12; % corresponds to 12*1.5 = 18km
    upper_height_index = 40; % corresponds to 40*1.5 = 60km
    
    tapering_length = 6; % corresponds to 6*1.5 = 9km window at edge
    
    Model.Tpg = sgwex_airs_apply_height_window(Model.Tpg,lower_height_index,upper_height_index,tapering_length);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% 8. remove altitude amplitude trend for the s-transform...
%     disp('- Removing altitude-amplitude trend...')
    
    ref_z = 40; % km
    H = scale_height(ref_z);
    sfvec = (exp(-(Model.Alt-ref_z) ./ (2*H)))';
    altampscaling = sfvec;
    sf = permute(repmat(sfvec,1,gsz(1),gsz(2)),[2 3 1]);
    Model.Tpg_sc = Model.Tpg .* sf;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% 9. 3DST!!!!
%     disp('- 3DST...')
    
    % first thing to do is fix NaNs. There will be loads from the AIRS
    % regridding.
    IN = Model.Tpg_sc;
    nanlocs = isnan(IN);
    IN(nanlocs) = 0;
    
    % Settings:
    nfreqs = 1000;
    minwavelengths = [25 25 3];
    c = [0.25 0.25 0.25];
    grid_spacing = [900 1200 75] ./ gsz;
    
% % % % %     % For the 3DST, apply a TEENY bit of horizontal smoothing to reduce noise in AIRS.
% % % % %     horz_fwhm = [13.5 13.5]; % km, airs footprint...
% % % % %     sig = (horz_fwhm ./ grid_spacing(1:2)) ./ 2.355;
% % % % %     IN = double(imgaussfilt(IN,sig,...
% % % % %         'padding','replicate',...
% % % % %         'filterdomain','spatial'));
    
    % time the 3DST...
    tic
    Model.ST = nph_ndst(IN,nfreqs,grid_spacing,c,'minwavelengths',minwavelengths);
    t1 = toc;
    disp(['3DST computed in ' num2str(round(t1)) 's.'])
    
    % remove altitude-amplitude scaling and replace nans...
    props = {'IN','HA','HR','A','R','BoostFactor','C'};
    for p = 1:length(props)
        Model.ST.(props{p}) = Model.ST.(props{p}) ./ sf;
        Model.ST.(props{p})(nanlocs) = NaN;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% 10. compute MF...
    
%     disp('- Computing MF...')
    
    N = 0.02; % gives BV period of about 5 minutes.
    g = 9.81; % m
    
    % create density and brunt vaisala frequency
    dens = permute(repmat(alt2dens(Model.Alt),gsz(1),1,gsz(2)),[1 3 2]); % kg/m^3
    
    % directional mf:
    Model.MF1 = 0.5 .* dens .* (g./N).^2 .* ((Model.ST.A)./Model.bgg).^2 .* (Model.ST.F1./Model.ST.F3);
    Model.MF2 = 0.5 .* dens .* (g./N).^2 .* ((Model.ST.A)./Model.bgg).^2 .* (Model.ST.F2./Model.ST.F3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% 11. Save!!!!
    
    % select properties you want:
    MA = struct;
    
    
    % make it look like AIRS:
    MA.T = single(Model.Ta_sm);
    MA.Lon = single(Airs.Lon);
    MA.Lat = single(Airs.Lat);
    MA.Alt = single(Model.Alt);
    MA.Time = Model.Time;
    
    % AIRS details...
    MA.AirsGranule = nearestgranule;
    MA.AirsFiles = Airs.Files;
    
    % gridded data:
    MA.Long = single(Model.Long);
    MA.Latg = single(Model.Latg);
    MA.Tg = single(Model.Tg);
    MA.Tpg = single(Model.Tpg);
    MA.bgg = single(Model.bgg);
    
    % s-transform:
    props = {'IN','numscales','scales','freqs','point_spacing','c','HA','HR','C','F1','F2','F3'};
    for p = 1:length(props)
        MA.ST.(props{p}) = Model.ST.(props{p});
    end
    
    % MF:
    MA.MF1 = Model.MF1;
    MA.MF2 = Model.MF2;
    
    MA.Name = 'Model as AIRS';
    
    Model = MA;
    
    
    % and save!
    disp(['- Saving to: ' savedirec granule '_Model_as_AIRS_3DST_AG.mat'])
    save([savedirec granule '_Model_as_AIRS_3DST_AG.mat'],'Model');
    disp('Saved!')
    
    
    
end


% end % END FUNCTION





















% 
% 
% return
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%% 1 ROTATE MODEL =========================================================
% 
% disp('rotating model...')
% 
% % True (rough) shape of model domain on square lat/lon coords:
% %        4--__    __--3
% %       /     ----     \
% %      /                \
% %     /                  \
% %    /                    \
% %   1-__                __-2
% %       --____________--
% 
% Model.Time = datestr(Model.MatlabTime);
% 
% % xy grid of the model in km (1200x900km)
% x = linspace(-600,598.5,800);
% y = linspace(-450,448.5,600);
% [D1,D2] = ndgrid(x,y);
% 
% % take rotate pole centre lat/lon
% rclat = double(Model.Lat(300));
% rclon = double(Model.Lon(400));
% 
% % rotate to true centre lat/lon
% clat = rclat - (90-Model.NorthPoleLatLon(1));
% clon = wrapTo180(rclon - (180-Model.NorthPoleLatLon(2)));
% 
% % find azimuth and distance to of every point from centre:
% az = atan2d(D1,D2);
% d = quadadd(D1,D2);
% 
% % use reckon to compute the lat/lon of every point on the grid:
% [Model.Lat,Model.Lon] = reckon(clat,clon,km2deg(d),az);
% 
% 
% % return
% 
% % new order?
% 
% % 1. fix rotate pole
% % 2. find unique vert_res's
% % 3. for each vert res, apply this to model, with footprint smoothing too.
% % 4. now you just need to evaluate model at airs position and height.
% 
% 
% 
% %% okay, now do the thing:
% % first, put the airs and model on the same vertical grid. There might be
% % some aliasing effects here but I don't think they'll be important. I
% % really hope I'm right!! considering the vertical scales of the waves in
% % the stratosphere I don't think a bit of aliasing will matter.
% % Airs.
% 
% % NOPE
% 
% % do the horizontal regridding first.
% 
% % NOPE
% 
% % NOPE
% 
% % do the vertical regridding first. You need it for applying the vertical
% % resolution.
% % if we simply regrid to 118 evenly spaced vertical levels, we get a
% % spacing of 650m from the surface to 75km. This is not bad. I don't expect
% % waves to change much in the vertical in the stratosphere over 650m,
% % (certainly not those large amplitudes visible to airs), so this is
% % probably reasonable. Hopefully aliasing shouldn't be an issue.
% 
% Model.T = permute(Model.temp,[2 1 3]);
% Model.Alt = Model.Alt_1 ./ 1000;
% % quickly regrid in the vertical before the vertical smoothing:
% sz = size(Model.T);
% F = griddedInterpolant({single(1:600),single(1:800),Model.Alt},Model.T,'linear','none');
% % compute vertical grid spacing:
% even_alt = double(linspace(min(Model.Alt),max(Model.Alt),sz(3)));
% vert_grid_spacing = nanmean(diff(even_alt));
% % evaluate the new gridded temp:
% Model.T_vert_grid = F({1:600,1:800,even_alt});
% 
% % get ready to simulate the airs footprint:
% horz_fwhm = [13.5 13.5]; % km
% horz_grid_spacing = [1.5 1.5]; % km
% 
% % great, now apply the AIRS vertical resolution. We have to do this now
% % because we can't cope with NaNs, which will inevitably arise drom the
% % regridding onto the AIRS scan track.
% %
% vert_resi = floor(airs3d_vert_res(even_alt,'nearest')); % give benefit of the doubt ;)
% uvr = unique(vert_resi);
% 
% Model.T_sm = nan(size(Model.T_vert_grid));
% 
% 
% for i = 1:length(uvr)
%     
%     disp(uvr(i))
%     
%     % compute 3d fwhm and grid spacing:
%     fwhm = [horz_fwhm uvr(i)];
%     grid_spacing = [horz_grid_spacing vert_grid_spacing];
%     
%     % smooth the model for this vert_res:
%     sig = (fwhm ./ grid_spacing) ./ 2.355;
%     smoo = imgaussfilt3(Model.T,sig,...
%         'padding','replicate',...
%         'filterdomain','spatial');
%     
%     % identify all the altitudes with this vert res:
%     zinds = find(vert_resi == uvr(i));
%     
%     % and assign the smoothed values:
%     Model.T_sm(:,:,zinds) = smoo(:,:,zinds);
%     
%     
% end
% 
% 
% return
% 
% 
% 
% 
% % % % %
% % % % %
% % % % % % new approach: do it for each AIRS height and regrid as you go. A
% % % % % % little slower but a lot simpler. And fairer I guess.
% % % % %
% % % % % % Model.T_sm = nan(size(Model.T_vert_grid));
% % % % %
% % % % % % just use the vert_res from each airs layer
% % % % % vert_res = airs3d_vert_res(Airs.Alt);
% % % % %
% % % % % % get ready to simulate the airs footprint:
% % % % % horz_fwhm = [13.5 13.5]; % km
% % % % % horz_grid_spacing = [1.5 1.5]; % km
% % % % %
% % % % % % create an interpolant for the horizontal model grid.
% % % % % F_model_horz = scatteredInterpolant(Model.Lon(:),Model.Lat(:),zeros(size(Model.Lon(:))),'linear','none');
% % % % %
% % % % % % and create a dummy interpolant for evaluating in the vertical:
% % % % % F_model_vert = griddedInterpolant({1:600,1:800,even_alt},zeros(size(Model.T)),'linear','none');
% % % % %
% % % % % Model.T_ag = nan(size(Airs.T));
% % % % %
% % % % % % now for each vertical layer:
% % % % % for z = 1:length(Airs.Alt)
% % % % %
% % % % %     disp(z)
% % % % %
% % % % %     % extract airs altitude for this layer
% % % % %     airs_z = Airs.Alt(z);
% % % % %
% % % % %     % deal with possible nans at the bottom layer
% % % % %     if airs_z == 0, airs_z = min(even_alt(:)); end
% % % % %
% % % % %     % compute 3d fwhm and grid spacing:
% % % % %     fwhm = [horz_fwhm vert_res(z)];
% % % % %     grid_spacing = [horz_grid_spacing vert_grid_spacing];
% % % % %
% % % % %     % smooth the model for this vert_res:
% % % % %     sig = (fwhm ./ grid_spacing) ./ 2.355;
% % % % %     smoo = imgaussfilt3(Model.T,sig,'padding','replicate');
% % % % %
% % % % %     % evaluate this smoothed temp at the appropriate height:
% % % % %     F_model_vert.Values = smoo;
% % % % %     T_layer = F_model_vert({1:600,1:800,airs_z});
% % % % %
% % % % %     % finally evaulate this layer onto the AIRS horz grid:
% % % % %     F_model_horz.Values = double(T_layer(:));
% % % % %     Model.T_ag(:,:,z) = F_model_horz(Airs.Lon,Airs.Lat);
% % % % %
% % % % % end
% % % % %
% % % % % % does it look as pretty? no. is it fairer? yes!
% % % % %
% % % % %
% 
% 
% 
% 
% return
% 
% 
% 
% 
% 
% 
% 
% % SGWEX MODEL CENTRE POINT:
% clat = -54.5269874182436;
% clon = -37.1367495437688;
% 
% 
% %% LOAD TOPO AND MAP ======================================================
% 
% if ~exist('tp','var')
%     disp('Loading Map and Topography...')
%     load('/Users/neil/Drive/MATLAB/topography/easy_tenth_degree_topography/Global.mat');
%     Ftp = griddedInterpolant({-90:0.1:90,-180:0.1:180},tp,'linear','none');
% end
% % if ~exist('Map','var')
% %     %     load('/Users/neil/Drive/MATLAB/imagery/SAAP_Region_Map.mat');
% % %     load('/Users/neil/Drive/MATLAB/imagery/world_map_light.mat');
% %     load('/Users/neil/Drive/MATLAB/imagery/map.mat');
% %     FMap = griddedInterpolant({Map.Lat,Map.Lon,1:3},single(Map.Map),'linear','none');
% % end
% 
% tpr = Ftp(MA.Model.Lat,MA.Model.Lon);
% tpr(tpr < 0 ) = 0;
% tpr = tpr ./ 1000;
% 
% 
% disp('Ready to plot.')
% 
% return
% 
% 
% 
% %%
% 
% %% ========================================================================
% %% PLOTTING!!!!!!!!! ======================================================
% %% ========================================================================
% 
% figure; hold all; grid on;
% whitefig;
% 
% % figpos([188         155        1376         921])
% figpos([0.75 0.7])
% 
% vert_gap = 0.06;       horz_gap = 0.065;
% lower_marg = 0.15;   upper_marg = 0.05;
% left_marg = 0.08;   right_marg = 0.08;
% 
% rows = 2; cols = 3;
% 
% subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);
% 
% 
% fs = 20;
% 
% sq = @squeeze;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % VERTICAL CROSS SECTIONS =================================================
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% latindex = 30;
% 
% a = sq(Airs.Tpg(latindex,:,:));
% m = sq(M.Model.Tp(:,latindex*10,:));
% % m = sq(M.Model.Tpg(latindex,:,:));
% ma = sq(MA.Model.Tpg(latindex,:,:));
% 
% [x,y] = ndgrid(linspace(-600,600,80),1.5:1.5:75);
% [xx,yy] = ndgrid(linspace(-600,600,800),M.Model.Alt./1000);
% 
% % Color levels:
% aclim = [-7 7];
% mclim = [-16 16];
% maclim = [-7 7];
% 
% a(a > max(aclim)) = max(aclim);
% a(a < min(aclim)) = min(aclim);
% m(m > max(mclim)) = max(mclim);
% m(m < min(mclim)) = min(mclim);
% ma(ma < min(maclim)) = min(maclim);
% ma(ma < min(maclim)) = min(maclim);
% 
% clev = 22;
% 
% % AIRS
% subplot(rows,cols,1)
% contourf(x,y,a,linspace(min(aclim),max(aclim),clev),'edgecolor','none');
% clim(aclim);
% ylabel('Altitude (km)')
% 
% % Alpha layers to smooth out the top and bottom tapering
% alphalevels = linspace(0,1,5);
% for i = 1:5
%     ztop = nan(80,12-i);
%     zbot = nan(80,15-i);
%     
%     hold on; ptop = pcolor(x(:,end-size(ztop,2)+1:end),y(:,end-size(ztop,2)+1:end),ztop); shat;
%     hold on; pbot = pcolor(x(:,1:size(zbot,2)),y(:,1:size(zbot,2)),zbot); shat;
%     
%     set(ptop,'facecolor','w','facealpha',alphalevels(i));
%     set(pbot,'facecolor','w','facealpha',alphalevels(i));
%     
% end
% 
% % MODEL
% subplot(rows,cols,2)
% contourf(xx,yy,m,linspace(min(mclim),max(mclim),clev),'edgecolor','none');
% clim(mclim);
% 
% % MODEL AS AIRS
% subplot(rows,cols,3)
% y(:,1) = 0;
% contourf(x,y,ma,linspace(min(maclim),max(maclim),clev),'edgecolor','none');
% clim(maclim);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % HORIZONTAL CROSS SECTIONS ===============================================
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% alt = 45;
% [~,zindex] = min(abs(Airs.Alt - alt));
% [~,mzindex] = min(abs(M.Model.Alt./1000 - alt));
% 
% a = sq(Airs.Tpg(:,:,zindex))';
% m = sq(M.Model.Tp(:,:,mzindex));
% % m = sq(M.Model.Tpg(:,:,zindex));
% ma = sq(MA.Model.Tpg(:,:,zindex))';
% 
% [x,y] = ndgrid(linspace(-600,600,80),linspace(-450,450,60));
% [xx,yy] = ndgrid(linspace(-600,600,800),linspace(-450,450,600));
% 
% a(a > max(aclim)) = max(aclim);
% a(a < min(aclim)) = min(aclim);
% m(m > max(mclim)) = max(mclim);
% m(m < min(mclim)) = min(mclim);
% ma(ma < min(maclim)) = min(maclim);
% ma(ma < min(maclim)) = min(maclim);
% 
% clev = 22;
% 
% % AIRS
% subplot(rows,cols,4)
% contourf(x,y,a,linspace(min(aclim),max(aclim),clev),'edgecolor','none');
% clim(aclim);
% ylabel('Meridional Distance (km)')
% 
% % MODEL
% subplot(rows,cols,5)
% contourf(xx,yy,m,linspace(min(mclim),max(mclim),clev),'edgecolor','none');
% clim(mclim);
% 
% % MODEL AS AIRS
% subplot(rows,cols,6)
% contourf(x,y,ma,linspace(min(maclim),max(maclim),clev),'edgecolor','none');
% clim(maclim);
% 
% % Formatting
% letters = alphabet(1,6);
% titles = {'AIRS' 'Model' 'Model as AIRS'};
% 
% % COASTLINE:
% C = nph_draw_coastline([min(MA.Model.Lon(:)) min(MA.Model.Lat(:)) ; max(MA.Model.Lon(:)) max(MA.Model.Lat(:))],0,'noplot');
% % convert to x y:
% for i = 1:length(C)
%     [arclen,az] = distance(clat,clon,C(i).Lat,C(i).Lon);
%     C(i).x = deg2km(arclen) .* sind(az);
%     C(i).y = deg2km(arclen) .* cosd(az);
% end
% 
% 
% % PLOT!!!!
% 
% for ax = 1:6
%     
%     subplot(rows,cols,ax); grid on;% box off;
%     setfont(fs)
%     
%     set(gca,'tickdir','in','linewi',1.25)
%     
%     % Colormap
%     cmap = colormap(cbrew('RdBu',64));
%     
%     switch ax
%         case {1,2,3}
%             
%             % Ticks
%             set(gca,'ytick',0:15:75)
%             set(gca,'xtick',-600:200:600)
%             
%             xlim([-600 600]);
%             ylim([0 75]);
%             
%             % SG Topograpy
%             hold on; fill([linspace(-600,600,800) -600],[max(tpr) 0],'k')
%             
%             % Panel Letter
%             xlm = get(gca,'xlim');
%             ylm = get(gca,'ylim');
%             text(xlm(1) + 0.035*abs(diff(xlm)),ylm(1) + 0.925*abs(diff(ylm)),['(' letters{ax} ')'],'fontsize',1.6*fs)
%             
%             % Title
%             t = title(titles{ax});
%             t.Position = t.Position .* [1 1.025 1];
%             
%             % Model sponge Layer
%             switch ax
%                 case {2,3}
%                     xl = get(gca,'xlim');
%                     hold on; plot(xl,[58.5 58.5],'color',rgbtrip(.5),'linest','--','linewi',1.5);
%                     %             if ax == 3
%                     text(1.035*xl(2),58.5,{'Model','Sponge','Layer'},'fontsize',0.65*fs,'color',rgbtrip(.35));
%                     %             end
%             end
%             
%             
%             % HORZ SLICES
%         case {4,5,6}
%             
%             xlim([-600 600]);
%             ylim([-450 450]);
%             
%             % Ticks
%             set(gca,'ytick',-400:200:400)
%             set(gca,'xtick',-600:200:600)
%             
%             xlabel('Zonal Distance (km)')
%             
%             % SG Outline:
%             for i = 1:length(C)
%                 hold on; plot(C(i).x,C(i).y,'k','linewi',1.25)
%             end
%             
%             % Colorbar
%             Cbars.(['cbar' num2str(ax)]) = colorbar('location','southoutside');
%             apos = get(gca,'position');
%             Cbars.(['cbar' num2str(ax)]).Position(3) = 0.75*apos(3);
%             Cbars.(['cbar' num2str(ax)]).Position(1) = apos(1) + (0.5*(apos(3) - Cbars.(['cbar' num2str(ax)]).Position(3))) - 0.001;
%             Cbars.(['cbar' num2str(ax)]).Position(2) = Cbars.(['cbar' num2str(ax)]).Position(2) - 0.13;
%             
%             % T' (K)
%             tpos = Cbars.(['cbar' num2str(ax)]).Position;
%             tpos(1) = tpos(1) + tpos(3) + 0.005;
%             %     tpos(2) = tpos(2) - 0.01;
%             annotation('textbox','position',tpos,'string','T'' (K)','fontsize',fs,'linest','none','verticalalignment','middle')
%             
%             % Panel Letter
%             xlm = get(gca,'xlim');
%             ylm = get(gca,'ylim');
%             tlet = text(xlm(1) + 0.035*abs(diff(xlm)),ylm(1) + 0.925*abs(diff(ylm)),['(' letters{ax} ')'],'fontsize',1.6*fs);
%             
%             % Altitude marker:
%             text(tlet.Position(1)+150,tlet.Position(2),['z = ' num2str(round(alt)) 'km'],'fontsize',fs,'horizontalalignment','left')
%             
%             
%     end
% end
% 
% Cbars.cbar4.Ticks = -8:2:8;
% Cbars.cbar5.Ticks = -16:4:16;
% Cbars.cbar6.Ticks = -8:2:8;
% 
% % fix AIRS clims:
% Cbars.cbar4.Limits = [-7 7];
% 
% 
% 
% return
% 
% %% EXPORT? ================================================================
% 
% savename = ['~/Desktop/' granule '_CrossSections'];
% 
% disp(['Exporting to "' savename '"...'])
% nph_saveas(gcf,savename,'png')
% disp('Done.')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
