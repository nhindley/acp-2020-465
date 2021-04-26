

%%%% EDIT: NEW APPROACH : SHOW AIRS AND MODEL AS AIRS ON THE SAME AIRS
%%%% GRID USING THE MOSAIC APPROACH TO SHOW THE AIR FOOTPRINTS



% 3x3 panel.
% Vertical and Horizontal corss sections through an AIRS granule, the SGWEX
% Model, and the Model as Airs.


%% DIRECTORIES ============================================================

matlabdirec = '/Users/neil/Drive/MATLAB/';
airsdirec = '/Users/neil/data/';
modeldirec = '/Users/neil/data/';
modelasairsdirec = '/Users/neil/data/';


%% LOAD AIRS AND MODEL ====================================================
% ---------------------------------------------------- single granules ----

% granule = '20150705-0300';
granule = '20150705-0300';

disp(['Loading AIRS and Model for ' granule '...'])

% load AIRS
load([airsdirec granule '_AirsSG.mat'])

% load full res model:
load([modeldirec granule '_Model_3DST.mat'])
% ^ try again with old model output. It's on the 1.5km vertical grid which
% is not ideal but never mind. second best solution today

disp('Loading Topography...')
load('/Users/neil/Drive/MATLAB/topography/SRTM/srtm_29_23/srtm_29_23.mat')

% trim to island area
topo = double(Topo.Topo(1:1500,2000:5500)) ./ 1000;
topolat = Topo.Lat(1:1500);
topolon = Topo.Lon(2000:5500);
% fix sea level
topo(topo < 0) = 0;

% take max:
topoline = nanmax(topo,[],1);

% get this is distance:
clon = Model.Lon(300,400);
clat = Model.Lat(300,400);
[d,az] = distance(clat,clon,clat,topolon);
topoxdist = deg2km(d).*sind(az);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MODEL -> MODEL AS AIRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. simulate the airs horz footrpint in the model. Nice and easy.
horz_fwhm = [13.5 13.5]; % km
horz_grid_spacing = [1.5 1.5]; % km
sig = (horz_fwhm ./ horz_grid_spacing) ./ 2.355;
Model.T_sm = double(imgaussfilt(Model.T,sig,...
    'padding','replicate',...
    'filterdomain','spatial'));

% 2. regrid to the AIRS scantrack. to do this, you need to get distance
% of all airs points from model centre, then you can use girdded interpolant for the model.

% 2.1 AIRS distance from model centre
clon = Model.Lon(300,400);
clat = Model.Lat(300,400);
[d,az] = distance(clat,clon,Airs.Lat,Airs.Lon);
xdist = deg2km(d).*sind(az);
ydist = deg2km(d).*cosd(az);

% 2.2 now you only need to make a *gridded* interpolant for the model:
mx = linspace(-600,598.5,800);
my = linspace(-450,448.5,600);
mz = Model.Alt;
Fh_model = griddedInterpolant({my,mx,mz},Model.T_sm,'linear','none');

% 2.3. regrid to the AIRS scantrack
disp('regridding model to airs...')
szm = size(Model.T);
sza = size(Airs.T);
Model.Ta = nan(sza(1),sza(2),szm(3));
% make some big versions of the grid:
X = repmat(xdist,1,1,szm(3));
Y = repmat(ydist,1,1,szm(3));
Z = permute(repmat(Model.Alt(:),1,sza(1),sza(2)),[2 3 1]);
Model.Ta = Fh_model(Y,X,Z);

% 3. apply the airs vertical resolution. this should be faster now because
% 90x270x50 is a lot smaller than 600x800x118!
disp('apply airs vertical resolution(s)...')
vert_resi = floor(0.5*airs3d_vert_res(Model.Alt,'nearest')); % give benefit of the doubt ;)
vert_grid_spacing = nanmean(diff(Model.Alt));
uvr = unique(vert_resi);
Model.Ta_sm = zeros(size(Model.Ta));

% fix nans resulting from the airs regridding:
nanlocs = isnan(Model.Ta);
for z = 1:szm(3)
    layer = Model.Ta(:,:,z);
    layer(isnan(layer)) = nanmedian(linearise(layer));
    Model.Ta(:,:,z) = layer;
end

% now for each vertical resolution:
for i = 1:length(uvr)
    
    disp(uvr(i))
    
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
    
    % identify all the altitudes with this vert res:
    zinds = find(vert_resi == uvr(i));
    
    % and assign the smoothed values:
    Model.Ta_sm(:,:,zinds) = smoo(:,:,zinds);
    
end

% put back the nans for consistency:
Model.Ta_sm(nanlocs) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 4. ADD AIRS NOISE TO THE MODEL AS AIRS
disp('adding AIRS noise to the model...')
% 20150620-1700 looks good.
% best to use a real granule cos the noise is different for each altitude.
A = load(['/Volumes/SDBlue/data/sgwex/AIRS/three_corner_overlapping_granule_pairs/matlab/20150620-1700_AirsSG.mat']);
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
% put the airs noise on the regular vertical grid:
FF = griddedInterpolant({1:90,1:270,A.Airs.Alt},airsnoise,'linear','none');
airsnoise = FF({1:90,1:270,Model.Alt});

% and add it to the model as airs:
Model.Ta_sm = Model.Ta_sm + (0.5*airsnoise);

% and you're done!



% done!



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ASSEMBLING PROPERTIES TO PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('extracting perturbations...')

% remove background:
[Model.Tp,Model.bg] = nph_airs_4dp_detrend(Model.T,2,4);
[Model.Tpa_sm,Model.bg_sm] = nph_airs_4dp_detrend(Model.Ta_sm,2,4);
[Airs.Tp,Airs.bg] = nph_airs_4dp_detrend(Airs.T,1,4);

% compute model distance grids:
[D1,D2] = meshgrid(mx,my);

disp('evaluating vertical cuts...')

% ALONG CROSS TRACK SCAN?
% % Now evaluate the vertical cuts:
% xcut = 115;
% 
% Airs.Tpv = sq(Airs.Tp(:,xcut,:))';
% Model.Tpav_sm = sq(Model.Tpa_sm(:,xcut,:))';
% 
% % model is a little more tricky...
% F = griddedInterpolant({my,mx,mz},double(Model.Tp),'linear','none');
% 
% ystripi = linspace(min(ydist(:,xcut)),max(ydist(:,xcut)),1200)';
% ex = repmat(xstripi,1,length(Model.Alt));
% ey = repmat(ystripi,1,length(Model.Alt));
% ez = repmat(Model.Alt,length(xstripi),1);
% Model.Tpv = F(ey,ex,ez)';

% OR ALONG A ZONAL LINE?
ycut = -17;
% for ycut = 0:30
xstripi = linspace(-600,600,800);
ystripi = ycut .* ones(size(xstripi));

% airs:
ynudge = 0;
Airs.Tpv = nan(length(Airs.Alt),length(xstripi));
Fa = scatteredInterpolant(ydist(:),xdist(:),zeros(size(xdist(:))),'linear','none');
for z = 1:length(Airs.Alt)
    Fa.Values = linearise(Airs.Tp(:,:,z));
    Airs.Tpv(z,:) = Fa({ycut+ynudge,xstripi});
end
% % or another way with airs? nearest footrpint approach?
% [I1,I2] = ndgrid(1:90,1:270);
% FI1 = scatteredInterpolant(ydist(:),xdist(:),I1(:),'nearest','none');
% FI2 = scatteredInterpolant(ydist(:),xdist(:),I2(:),'nearest','none');
% i1 = FI1({ycut,xstripi});
% i2 = FI2({ycut,xstripi});
% % now evaluate along this line:
% linds = sub2ind(size(Airs.Tp(:,:,16)),i1,i2);
% for z = 1:length(Airs.Alt)
%     layer = Airs.Tp(:,:,z);
%     Airs.Tpv(z,:) = layer(linds);
% end

% model:
ynudge = 0;
Model.Tpv = nan(length(Model.Alt),length(xstripi));
Fm = griddedInterpolant({my,mx},zeros(600,800),'linear','none');
for z = 1:length(Model.Alt)
    Fm.Values = Model.Tp(:,:,z);
    Model.Tpv(z,:) = Fm({ycut+ynudge,xstripi});
end

% model as airs:
ynudge = 0;
Model.Tpav_sm = nan(length(Model.Alt),length(xstripi));
for z = 1:length(Model.Alt)
    Fa.Values = linearise(Model.Tpa_sm(:,:,z));
    Model.Tpav_sm(z,:) = Fa({ycut+ynudge,xstripi});
end

% % squish em good:
% Model.Tpv       = sq(nanmedian(Model.Tpv,2));
% Airs.Tpv        = sq(nanmedian(Airs.Tpv,2));
% Model.Tpav_sm   = sq(nanmedian(Model.Tpav_sm,2));

% % and squish the strip:
% ystripi         = nanmean(ystripi) .* ones(size(xstripi));

% make some grids for us to plot these cuts with:
[x_model_as_airs,z_model_as_airs] = meshgrid(xstripi,Model.Alt);
[x_model,z_model] = meshgrid(xstripi,Model.Alt);

% finally, stick the AIRS cut on a regular height grid like the model.
% It's just easier to look at, and it's only interpolation.
Fa = griddedInterpolant({Airs.Alt,xstripi},Airs.Tpv,'linear','none');
Airs.Tpvm = Fa({Model.Alt,xstripi});
[x_airs,  z_airs] = meshgrid(xstripi,Model.Alt);

% cla
% hold on; contourf(xstripi,Model.Alt,Model.Tpv,31,'edgecolor','none'); shat; clim([-10 10]); colormap(cbrew)
% title(num2str(ycut))
% drawnow;
% 
% 
% return
% 
% figure; pcolor(Airs.Tpv); shat; clim([-10 10]); colormap(cbrew)
% figure; pcolor(Model.Tpv); shat; clim([-15 15]); colormap(cbrew)
% figure; pcolor(Model.Tpav_sm); shat; clim([-10 10]); colormap(cbrew)


disp('creating topography profile...')




return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PLOTTY WOTTY PLOT PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([1 0.45])

%-------------------------------------------------------
vert_gap    = 0.05;    horz_gap    = 0.06;
lower_marg  = 0.25;    upper_marg  = 0.05;
left_marg   = 0.05;    right_marg  = 0.05;

rows = 1; cols = 3;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 22;

%--------------------------------------------------------

part = 3;

drawmosaic = 0;

airszlev = 16; % 45km
modelzlev = 30; % 45km

% define boundaries
% xbox = [-42.1 -29.9];
% ybox = [-58.1 -50.5];


%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch part
    case 1
        markersizes = {...
            4.5,...
            NaN,...
            4.5}; % PART 1
        
        box = {...
            (-37.5)+plusminus(13),(-54.5)+plusminus(6.5);...
            (-37.5)+plusminus(12),(-54.5)+plusminus(5.5);...
            (-37.5)+plusminus(12),(-54.5)+plusminus(5.5)}; % PART 1
        
        TBP = {...
            Airs.T(:,:,airszlev)-9,...
            Model.T(:,:,modelzlev),...
            Model.Ta_sm(:,:,modelzlev)}; % PART 1
        
        tbpxy = {...
            Airs.Lon,Airs.Lat;...
            Model.Lon,Model.Lat;...
            Airs.Lon,Airs.Lat}; % PART 1
        
        mosaics = {
            1,...
            0,...
            1}; % PART 1
        
        clims = {...
            235+plusminus(12.5),...
            235+plusminus(12.5),...
            235+plusminus(12.5)}; % PART 1
        
        cbarlabels = {...
            'Temperature (K)',...
            'Temperature (K)',...
            'Temperature (K)'}; % PART 1
        
        cbarticks = {...
            200:5:300,...
            200:5:300,...
            200:5:300}; % PART 1
        
        cbarlims = {...
            clims{1},...
            clims{2},...
            clims{3}}; % PART 1
        
        %%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
        markersizes = {...
            7.5,...
            NaN,...
            7.5}; % PART 2
        
        box = {...
            [-500 500],plusminus(400);...
            plusminus(500),plusminus(400);...
            plusminus(500),plusminus(400)}; % PART 2
        
        TBP = {...
            0.9*Airs.Tp(:,:,airszlev),...
            Model.Tp(:,:,modelzlev),...
            Model.Tpa_sm(:,:,modelzlev)}; % PART 2
        
        tbpxy = {...
            xdist,ydist;...
            D1,D2;...
            xdist,ydist}; % PART 2
        
        mosaics = {
            1,...
            0,...
            1}; % PART 2
        
        clims = {...
            plusminus(8),...
            plusminus(12),...
            plusminus(8)}; % PART 2
        
        cbarlabels = {...
            'GW Perturbations (K)',...
            'GW Perturbations (K)',...
            'GW Perturbations (K)'}; % PART 2
        
        cbarticks = {...
            -15:5:15,...
            -15:5:15,...
            -15:5:15}; % PART 2
        
        cbarlims = {...
            clims{1},...
            clims{2},...
            clims{3}}; % PART 2
        
    %%%%% PART 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
        markersizes = {...
            8,...
            NaN,...
            8}; % PART 3
        
        box = {...
            plusminus(500),[0 70];...
            plusminus(500),[0 70];...
            plusminus(500),[0 70]}; % PART 3
        
        TBP = {...
            Airs.Tpvm,...
            Model.Tpv,...
            Model.Tpav_sm}; % PART 3
        
        tbpxy = {...
            x_airs,z_airs;...
            x_model,z_model;...
            x_model_as_airs,z_model_as_airs}; % PART 3
        
        mosaics = {
            0,...
            0,...
            0}; % PART 3
        
        clims = {...
            plusminus(8),...
            plusminus(12),...
            plusminus(8)}; % PART 3
        
        cbarlabels = {...
            'GW Perturbations (K)',...
            'GW Perturbations (K)',...
            'GW Perturbations (K)'}; % PART 3
        
        cbarticks = {...
            -20:2:20,...
            -20:4:20,...
            -20:2:20}; % PART 3
        
        cbarlims = {...
            clims{1},...
            clims{2},...
            clims{3}}; % PART 3
        
        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PLOTTY PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ax = 1:3
    
    subplot(rows,cols,ax)
    
    axx = gca;
    
    xbox = box{ax,1};
    ybox = box{ax,2};
    
    % set limits first to speed up mosaic
    xlim(xbox)
    ylim(ybox)
    
    % set model alt to go from zero (you muppet)
    if part == 3
        tbpxy{ax,2}(1,:) = 0;
        tbpxy{ax,2}(1,:) = 0;
        tbpxy{ax,2}(1,:) = 0;
    end
    
    % and set tick marks first:
    switch part
        case 1
            xtick(-80:4:0);
            xminortick(-80:2:0)
            ytick(-80:2:-30);
            yminortick(-80:1:-30);
        case 2
            xtick(-600:200:600);
            xminortick(-600:50:600)
            ytick(-600:100:600);
            yminortick(-600:50:600);
        case 3
            xtick(-600:200:600);
            xminortick(-600:50:600)
            ytick(0:10:70);
            yminortick(0:5:70);
    end
    
    % some manual grid lines are needed first:
    gridlinespec = {'linewi',1.5,'color',[0.15 0.15 0.15 0.15]};
    gl.xx = repmat(axx.XTick,2,1);
    gl.xy = repmat(axx.YLim,length(axx.XTick),1)';
    gl.yy = repmat(axx.YTick,2,1);
    gl.yx = repmat(axx.XLim,length(axx.YTick),1)';
    hold on; plot(gl.xx,gl.xy,gridlinespec{:});
    hold on; plot(gl.yx,gl.yy,gridlinespec{:});
    
    %%%% DATA!!!
    if drawmosaic
        switch mosaics{ax}
            case 1
                switch part
                    case 1
                        switch ax
                            case 1
                                % grey patch underneath:
                                hold on; P = patch(...
                                    Airs.Lon(edge_indeces(Airs.Lon)),...
                                    Airs.Lat(edge_indeces(Airs.Lat)),...
                                    rgbtrip(.975));
                                P.EdgeColor = 'none';
                            case 3
                                % grey patch underneath:
                                hold on; P = patch(...
                                    Model.Lon(edge_indeces(Model.Lon)),...
                                    Model.Lat(edge_indeces(Model.Lat)),...
                                    rgbtrip(.975));
                                P.EdgeColor = 'none';
                        end
                    case 2
                        % grey patch underneath:
                        set(gca,'color',rgbtrip(0.975))
                end
                % MOSAIC PLOT
                markersize = markersizes{ax};
                data = TBP{ax};
                cl = clims{ax}; clim(cl); nclev = 31;
                cmap = cbrew('RdBu',nclev);
                Fcmap = griddedInterpolant({linspace(cl(1),cl(2),nclev),1:3},cmap,'linear','nearest');
                for j = 1:size(data,2)
                    for i = 1:size(data,1)
                        % if in the region:
                        if inrange(tbpxy{ax,1}(i,j),xbox) && inrange(tbpxy{ax,2}(i,j),ybox)
                            if ~isnan(data(i,j))
                                % get dot color:
                                dotcolor = Fcmap({data(i,j),1:3});
                                % draw a dot:
                                hold on; plot(tbpxy{ax,1}(i,j),tbpxy{ax,2}(i,j),'marker','o','markerfacecolor',dotcolor,'markersize',markersize,'markeredgecolor','none');
                            end
                        end
                    end
                    drawnow;
                end
            otherwise
                % STANDARD PCOLOR
                hold on; pcolor(tbpxy{ax,1},tbpxy{ax,2},TBP{ax}); shat;
        end
    else
        switch part
            case {1,2}
                
                % STANDARD PCOLOR
                hold on; pcolor(tbpxy{ax,1},tbpxy{ax,2},TBP{ax}); shat;
                
            case 3
                
                % STANDARD CONTOURF
                data = TBP{ax};
                data = imgaussfilt(data,[0.05 1]);
                % Alpha layers to smooth out the top and bottom tapering
%                 fading_rate = 0.9;
                switch ax
                    case 1 % AIRS
                        fading_rate = 0.9;
                        toplayers = find(inrange(Model.Alt,[58 100]));
                        for z = toplayers
                            fading_rate = fading_rate.^2;
                            data(z:end,:) = fading_rate.*data(z:end,:);
                        end
                        fading_rate = 0.9;
                        bottomlayers = find(inrange(Model.Alt,[0 20]));
                        for z = reverse(bottomlayers)
                            fading_rate = fading_rate.^2;
                            data(1:z,:) = fading_rate.*data(1:z,:);
                        end
                        % set final limits:
                        data(~inrange(Model.Alt,[15 62]),:) = 0;
                        %                     set(gca,'color',rgbtrip(.975))
                    case 3 % MODEL AS AIRS
%                         fading_rate = 0.9;
%                         toplayers = find(inrange(Model.Alt,[58 100]));
%                         for z = toplayers
%                             fading_rate = fading_rate.^2;
%                             data(z:end,:) = fading_rate.*data(z:end,:);
%                         end
                        fading_rate = 0.95;
                        bottomlayers = find(inrange(Model.Alt,[0 15]));
                        for z = reverse(bottomlayers)
                            fading_rate = fading_rate.^2;
                            data(1:z,:) = fading_rate.*data(1:z,:);
                        end
                        
                end
                
                nclev = 31;
                cl = clims{ax};
                clev = linspace(cl(1),cl(2),nclev);
                data = nph_fix2clims(data,cl);
                hold on; contourf(tbpxy{ax,1},tbpxy{ax,2},data,clev,'edgecolor','none');
        end
    end
    
    %%%% FORMATTING
    if part == 1
        axx.XTickLabel = lonlabels(axx.XTick);
        axx.YTickLabel = latlabels(axx.YTick);
    end
    
    % COLORBAR
    if part ~= 2
        switch ax
            case {1,2,3}
                
                cbar = colorbar('southoutside');
                cbar.Label.String = cbarlabels{ax};
                cbar.Label.Rotation = 0;
                %             cbar.Label.VerticalAlignment = 'middle';
                %             cbar.Label.HorizontalAlignment = 'left';
                cbar.TickDirection = 'out';
                cbar.Ticks = cbarticks{ax};
                cbar.LineWidth = 1.5;
                cbar.Limits = cbarlims{ax};
                
                % minortick
                cbar.Ruler.MinorTick = 'on';
                switch part
                    case 1
                        cbar.Ruler.MinorTickValues = 200:5:300;
                    case 3
                        cbar.Ruler.MinorTickValues = -20:20;
                end
                
                % reposition
                cbar.Position = [1 1 1 1] .* cbar.Position;
                drawnow
                apos = get(gca,'position');
                cbar.Position(3) = 0.6*apos(3);
                cbar.Position(1) = apos(1) + 0.2*apos(3);
                cbar.Position(2) = cbar.Position(2) - 0.165;
                cbar.Position(4) = 1.5*cbar.Position(4);
        end
    end
    
    % COASTLINE:
    switch part
        case 1
            C = nph_draw_coastline([xbox(1) ybox(1) ; xbox(2) ybox(2)],0,'noplot');
            for d = 1:length(C)
                if length(C(d).Lon) > 1
                    hold on; plot(C(d).Lon,C(d).Lat,'w','linewi',3)
                    hold on; plot(C(d).Lon,C(d).Lat,'k','linewi',2)
                end
            end
        case 2
            C = nph_draw_coastline([min(Model.Lon(:)) min(Model.Lat(:)) ; max(Model.Lon(:)) max(Model.Lat(:))],0,'noplot');
            for d = 1:length(C)
                if length(C(d).Lon) > 1
                    [d,az] = distance(clat,clon,C(d).Lat,C(d).Lon);
                    xd = deg2km(d).*sind(az);
                    yd = deg2km(d).*cosd(az);
                    hold on; plot(xd,yd,'w','linewi',3)
                    hold on; plot(xd,yd,'k','linewi',2)
                end
            end
        case 3
            % topography profile:
            hold on; P = fill(topoxdist,topoline,'k');
    end
    
    % draw model box:
    if part == 1
        linespec = {'color',mcolor(1),'linest','--','linewi',2};
        hold on; plot(Model.Lon(edge_indeces(Model.Lon)),Model.Lat(edge_indeces(Model.Lat)),linespec{:},'color','w','linest','-');
        hold on; plot(Model.Lon(edge_indeces(Model.Lon)),Model.Lat(edge_indeces(Model.Lat)),linespec{:});
    end
    
    % and some cross section lines
    switch part
        case 2
            linespec = {'linewi',2.5,'linest','--','color',[1 0 1 0.8]};
            hold on; p = plot(xstripi,ystripi-8,linespec{:});
    end
    
    clim(clims{ax})
    colormap(gca,cbrew('RdBu',31))
    
    % titles and letters
    drawnow;
    nph_text([0.025 0.87],['(' alphabet(ax+((part-1)*3)) ')'],'fontsize',1.8*fs,'HorizontalAlignment','center','textborder','w');
    %     nph_text([0.1 0.865],tits{ax},'fontsize',1.6*fs,'fontweight','bold','HorizontalAlignment','center','textborder','w');
    
    %%%% LABELS
    switch part
        case 1
            
        case 2
            ylabel('Meridional Distance (km)');
            xlabel('Zonal Distance (km)');
        case 3
            ylabel('Altitude (km)');
    end
    
    % Model sponge Layer
    if part == 3
        switch ax
            case {2,3}
            hold on; plot(axx.XLim,[58.5 58.5],'color',rgbtrip(.5),'linest','--','linewi',2);
            text(1.035*axx.XLim(2),64,{'Model','Sponge','Layer'},'fontsize',0.75*fs,'color',rgbtrip(.15));
        end
    end
    
    setfont(fs)
    set(gca,'linewi',1.5,'layer','top','tickdir','out')
    
    
    % outer line
    hold on; plot(xbox([1 2 2 1 1]),ybox([1 1 2 2 1]),'color','k','linewi',1.5);
    
end % next axes



return




%% EXPORT? ================================================================

savename = ['~/Desktop/' granule '_CrossSectionsAirsGrid_pt' num2str(part)];

disp(['Exporting to "' savename '"...'])
nph_saveas(gcf,savename,'png')
disp('Done.')




return



































