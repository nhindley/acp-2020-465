

% Eight 3D plots of AIRS/Model as AIRS propeties...

%%%% AG VERSION

% Does it in two parts, first showing T' |T'| LH LZ, then showing MFX MFY W
% and CpH.

% -------------------------------------------------------------------------

% As a demonstration of the 3DST working nicely, plot a 3D plot of one of
% the SGWEX granules over south georgia in 3D. With my new found
% proficiency in isosurfaces, should look quite nice!

% Compute the 3DST again in-house to use the most up-to-date version. I've
% added little fixes like HA and HR since I did the SG runs.


granule = '20150705-1700';

%% LOAD AIRS AND MODEL AS AIRS ============================================

disp('Loading AIRS and model-as-AIRS...')

% load(['/Users/neil/Drive/MATLAB/SGWEX/' granule '_AirsSG_3DST_4dp.mat'])
% load(['/Users/neil/Drive/MATLAB/SGWEX/' granule '_Model_as_AIRS_3DST_4dp.mat'])
% load(['/Volumes/SDBlue/data/sgwex/3DST/Airs/' granule '_AirsSG_3DST_4dp.mat'])
% load(['/Volumes/SDBlue/data/sgwex/3DST/Model_as_AIRS/' granule '_Model_as_AIRS_3DST_4dp.mat'])

load(['/Volumes/SDBlue/data/sgwex/3DST/Airs/AG/' granule '_Airs_3DST_AG.mat'])
load(['/Volumes/SDBlue/data/sgwex/3DST/Model_as_AIRS/AG/' granule '_Model_as_AIRS_3DST_AG.mat'])


%% LOAD TOPO AND MAP ======================================================

if ~exist('tp','var')
    disp('Loading Map and Topography...')
    load('/Users/neil/Drive/MATLAB/topography/easy_tenth_degree_topography/Global.mat');
    Ftp = griddedInterpolant({-90:0.1:90,-180:0.1:180},tp,'linear','none');
end
if ~exist('Map','var')
    load('/Users/neil/Drive/MATLAB/imagery/SAAP_Region_Map.mat');
    FMap1 = griddedInterpolant({Map.Lat,Map.Lon},single(Map.Map(:,:,1)),'linear','none');
    FMap2 = griddedInterpolant({Map.Lat,Map.Lon},single(Map.Map(:,:,2)),'linear','none');
    FMap3 = griddedInterpolant({Map.Lat,Map.Lon},single(Map.Map(:,:,3)),'linear','none');
    
end

% SGWEX MODEL CENTRE POINT:
clat = -54.5269874182436;
clon = -37.1367495437688;

xbox = plusminus(700);
ybox = plusminus(550);

xvec = xbox(1):1.5:xbox(2);
yvec = ybox(1):2.5:ybox(2);

[Y,X] = ndgrid(yvec,xvec);

% Convert distance box to Lat Lon:
[LAT,LON] = reckon(clat,clon,km2deg(quadadd(X,Y)),atan2d(X,Y));

% pull down to region:
tpr = Ftp(LAT,LON) ./ 1000; tpr(tpr < 0) = 0;
Mapr = uint8(cat(3,FMap1(LAT,LON),FMap2(LAT,LON),FMap3(LAT,LON))) * 1.5; % brighten map

% Grayscale?
% Mapr = rgb2gray(Mapr);

% Coastline?
boundingbox = [min(LON(:)) min(LAT(:)) ; max(LON(:)) max(LAT(:))];
C = nph_draw_coastline(boundingbox,0,'noplot');

% smooth topography?
tpr = smoothn(tpr,[3 3]);


%% RUN 3DST AGAIN =========================================================

runagain = 1;

if runagain
    
    disp('Re-doing 3DSTs...')
    
    sz = size(Airs.Tpg);
    
    sampling_intervals = [900 1200 75] ./ sz;
    
    
    % remove exp trend with altitude:
    ref_z = 40; % km
    H = 7; % km
    sfvec = (exp(-(Airs.Alt-ref_z) ./ (2*H)))';
    altampscaling = sfvec;
    sf = permute(repmat(sfvec,1,sz(1),sz(2)),[2 3 1]);
    Airs.Tpg_sc  = Airs.Tpg  .* sf;
    Model.Tpg_sc = Model.Tpg .* sf;
    
    % Set NaNs from the extrapolation to be zero:
    airsnanlocs = isnan(Airs.Tpg_sc);
    modelnanlocs = isnan(Model.Tpg_sc);
    Airs.Tpg_sc(airsnanlocs) = 0;
    Model.Tpg_sc(modelnanlocs) = 0;
    
    % Set nfreqs for guided fourier mode:
    nfreqs = 1000;
    minwavelengths = [30 30 3];
    c = [0.25 0.25 0.25];
    
    tic
    Airs.ST = nph_ndst(Airs.Tpg_sc,nfreqs,sampling_intervals,c,'minwavelengths',minwavelengths);
    t1 = toc;
    disp(['AIRS 3DST computed in ' num2str(round(t1,2)) 's.'])
    
    tic
    Model.ST = nph_ndst(Model.Tpg_sc,nfreqs,sampling_intervals,c,'minwavelengths',minwavelengths);
    t1 = toc;
    disp(['AIRS 3DST computed in ' num2str(round(t1,2)) 's.'])
    
    % replace NaNs from interpolating
    Airs.Tpg_sc(airsnanlocs) = NaN;
    Model.Tpg_sc(modelnanlocs) = NaN;
    
    % remove altitude-amplitude scaling and replace nans...
    props = {'HA','HR','A','R','BoostFactor','C'};
    for p = 1:length(props)
        Airs.ST.(props{p}) = Airs.ST.(props{p}) ./ sf;
        Airs.ST.(props{p})(airsnanlocs) = NaN;
        Model.ST.(props{p}) = Model.ST.(props{p}) ./ sf;
        Model.ST.(props{p})(modelnanlocs) = NaN;
    end
    
%     % Save original
%     Airs.ST.IN_orig = Airs.ST.IN;
%     Model.ST.IN_orig = Model.ST.IN;
    
    disp('Done!')
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE MF AGAIN

disp('Computing MF...')

Airs.bgg = Airs.Tg - Airs.Tpg;
Model.bgg = Model.Tg - Model.Tpg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ASSUME PROPAGATION INTO THE WIND FOR MOUNTAIN WAVES

posinds = Airs.ST.F2 > 0;
Airs.ST.F2               = -abs(Airs.ST.F2);
Airs.ST.F1(posinds)      = -Airs.ST.F1(posinds);
Airs.ST.F3(posinds)      = -Airs.ST.F3(posinds);

posinds = Model.ST.F2 > 0;
Model.ST.F2               = -abs(Model.ST.F2);
Model.ST.F1(posinds)      = -Model.ST.F1(posinds);
Model.ST.F3(posinds)      = -Model.ST.F3(posinds);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     % First, avoid model wraparound issues from the 3DST:
%     wvec = 1 - linearise(exp(-(((Model.Alt)-75).^2) ./ (2*(20./2.355))));
%
%     Model.ST.IN     = Model.ST.IN .* permute(repmat(wvec,1,sz(1),sz(2)),[2 3 1]);
%     Model.ST.A      = Model.ST.A  .* permute(repmat(wvec,1,sz(1),sz(2)),[2 3 1]);
%     Model.ST.HA     = Model.ST.HA .* permute(repmat(wvec,1,sz(1),sz(2)),[2 3 1]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create density and brunt vaisala frequency matrices (use SI units, units
% of temp and wavenumber cancel out so don't worry about those)
N = 0.02; % gives BV period of about 5 minutes.
g = 9.81; % m
H = 7; % km
dens = 1.225 .* exp(-Airs.Alt./H);
Dens = permute(repmat(alt2dens(Airs.Alt),sz(1),1,sz(2)),[1 3 2]); % kg/m^3

% Airs
Airs.MF1 = 0.5 .* Dens .* (g./N).^2 .* (Airs.ST.A./Airs.bgg).^2 .* (Airs.ST.F1./Airs.ST.F3);
Airs.MF2 = 0.5 .* Dens .* (g./N).^2 .* (Airs.ST.A./Airs.bgg).^2 .* (Airs.ST.F2./Airs.ST.F3);

% Model
Model.MF1 = 0.5 .* Dens .* (g./N).^2 .* (Model.ST.A./Model.bgg).^2 .* (Model.ST.F1./Model.ST.F3);
Model.MF2 = 0.5 .* Dens .* (g./N).^2 .* (Model.ST.A./Model.bgg).^2 .* (Model.ST.F2./Model.ST.F3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CONVERT TO mPa

Airs.MF1 = 1000 * Airs.MF1;
Airs.MF2 = 1000 * Airs.MF2;

Model.MF1 = 1000 * Model.MF1;
Model.MF2 = 1000 * Model.MF2;



return



%% PLOTTING ===============================================================

%-------------------------------------------------------
figure; hold all; whitefig; figpos([0.9 1])

vert_gap   = 0.06;      horz_gap    = 0.01;
lower_marg = 0.14;      upper_marg  = 0.06;
left_marg  = 0.075;     right_marg  = 0.01;

rows = 2; cols = 4;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

fs = 18;
%--------------------------------------------------------

%%% First a regular distance grid:
vec1 = linspace(-450,450,90);
vec2 = linspace(-800,800,120);
vec3 = linspace(1.5,75,50);
[D1,D2,D3] = ndgrid(vec2,vec1,vec3);

% Now force propagataion into the wind ====================================
poslocs = Airs.ST.F2 > 0;
Airs.ST.F1(poslocs) = -Airs.ST.F1(poslocs);
Airs.ST.F2(poslocs) = -Airs.ST.F2(poslocs);
% Airs.ST.F3 = abs(Airs.ST.F3);

poslocs = Model.ST.F2 > 0;
Model.ST.F1(poslocs) = -Model.ST.F1(poslocs);
Model.ST.F2(poslocs) = -Model.ST.F2(poslocs);
% Model.ST.F3 = abs(Model.ST.F3);

% N = 0.02; % rad/s
% g = 9.81; % m/s
% H = 7; % km
% Dens = 1.225 .* exp(-D3./H); % D3 is in km, so this is fine.

% background:
Airs.bgg = Airs.Tg - Airs.Tpg;
Model.bgg = Model.Tg - Model.Tpg;

% Recalculate MF based on this assumption:
Airs.MF1  = 0.5 .* Dens .* (g./N).^2 .* (Airs.ST.A./Airs.bgg).^2  .* (Airs.ST.F1./abs(Airs.ST.F3));
Airs.MF2  = 0.5 .* Dens .* (g./N).^2 .* (Airs.ST.A./Airs.bgg).^2  .* (Airs.ST.F2./abs(Airs.ST.F3));
Model.MF1 = 0.5 .* Dens .* (g./N).^2 .* (Model.ST.A./Model.bgg).^2 .* (Model.ST.F1./abs(Model.ST.F3));
Model.MF2 = 0.5 .* Dens .* (g./N).^2 .* (Model.ST.A./Model.bgg).^2 .* (Model.ST.F2./abs(Model.ST.F3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CONVERT TO mPa

Airs.MF1 = 1000 * Airs.MF1;
Airs.MF2 = 1000 * Airs.MF2;

Model.MF1 = 1000 * Model.MF1;
Model.MF2 = 1000 * Model.MF2;

% =========================================================================

% compute ST wavelength, theta and phi:
Airs.ST.L       = 1./quadadd(Airs.ST.F1,Airs.ST.F2,Airs.ST.F3);
Airs.ST.LH      = 1./quadadd(Airs.ST.F1,Airs.ST.F2);
% Airs.ST.Theta   = atan2d(Airs.ST.F2,Airs.ST.F1);
% Airs.ST.Phi     = atan2d(Airs.ST.F3,quadadd(Airs.ST.F1,Airs.ST.F2));

Model.ST.L      = 1./quadadd(Model.ST.F1,Model.ST.F2,Model.ST.F3);
Model.ST.LH     = 1./quadadd(Model.ST.F1,Model.ST.F2);
% Model.ST.Theta  = atan2d(Model.ST.F2,Model.ST.F1);
% Model.ST.Phi    = atan2d(Model.ST.F3,quadadd(Model.ST.F1,Model.ST.F2));

% GIVE AIRS NOISE A LITTLE SMOOTH AT LOW ALTITUDES FIRST
altlim = 25;
A.IN(:,:,Airs.Alt < altlim) = smoothn(Airs.ST.IN(:,:,Airs.Alt < altlim),[5 5 1]);
M.IN(:,:,Airs.Alt < altlim) = smoothn(Model.ST.IN(:,:,Airs.Alt < altlim),[3 3 1]);

% A.A(:,:,Airs.Alt < altlim)  = smoothn(A.A(:,:,Airs.Alt < altlim),[5 5 1]);

% try some variations:
M.IN    = smoothn(Model.ST.IN,[3 3 1]);
M.A     = smoothn(Model.ST.A,[3 3 1]);
M.HA    = smoothn(Model.ST.HA,[3 3 1]);
M.HR    = smoothn(Model.ST.HR,[3 3 1]);
M.MF    = smoothn(quadadd(Model.MF1,Model.MF2),[3 3 1]);
M.MF1   = smoothn(Model.MF1,[3 3 1]);
M.MF2   = smoothn(Model.MF2,[3 3 1]);
M.LH    = smoothn(Model.ST.LH,[3 3 1]);

A.IN    = smoothn(A.IN,[3 3 1]);
A.A     = smoothn(Airs.ST.A,[3 3 1]);
A.HA    = smoothn(Airs.ST.HA,[3 3 1]);
A.HR    = smoothn(Airs.ST.HR,[3 3 1]);
A.MF    = smoothn(quadadd(Airs.MF1,Airs.MF2),[3 3 1]);
A.MF1   = smoothn(Airs.MF1,[3 3 1]);
A.MF2   = smoothn(Airs.MF2,[3 3 1]);
A.LH    = smoothn(Airs.ST.LH,[3 3 1]);

% % try the log of MF?
% neginds1 = A.MF1 < 0;
% A.MF1 = abs(log10(1000*A.MF1)); % log mPa
% A.MF1(neginds1) = -A.MF1(neginds1);

% TRIM WRAPAROUND
M.IN(:,1:15,25:50) = 0.5*M.IN(:,1:15,25:50);
% M.A(:,1:10,:) = NaN;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DEFINE PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

part = 1;
overlay = 1;

% CHOOSE PART
switch part
    
    %%%% PART 1
    case 1
        
        % define isosurfaces:
        iso_cvals = [2 5];
        iso_avals = [0.35 1];
        
        airs_iso_cvals      = [2 5];
        model_iso_cvals     = [1.25 4]; % part 1
        
        tbp_isosurfaces = {...
            A.IN,...
            A.HA .* sf,...
            A.HA .* sf,...
            A.HA .* sf,...
            M.IN,...
            M.HA .* sf,...
            M.HA .* sf,...
            M.HA .* sf}; % part 1
        
        tbp_cvals = {...
            [-2.25 2.25],...
            airs_iso_cvals,...
            airs_iso_cvals,...
            airs_iso_cvals,...
            [-1.6 1.6],...
            model_iso_cvals,...
            model_iso_cvals,...
            model_iso_cvals}; % part 1
        
        tbp_avals = {...
            [1 1],...
            iso_avals,...
            iso_avals,...
            iso_avals,...
            [1 1],...
            iso_avals,...
            iso_avals,...
            iso_avals}; % part 1
        
        tbp_data = {...
            A.IN ./ sf,...
            A.A,...
            Airs.ST.LH,...
            abs(1./Airs.ST.F3),...
            M.IN ./ sf,...
            M.A,...
            Model.ST.LH,...
            abs(1./Model.ST.F3)}; % part 1
        
        cmaps = {...
            cbrew('RdBu',64),...
            cbrew('OrRd',64),...
            flipud(cbrew('Blues',64)),...
            cbrew('Greens',64),...
            cbrew('RdBu',64),...
            cbrew('OrRd',64),...
            flipud(cbrew('Blues',64)),...
            cbrew('Greens',64)}; % part 1
        
        clims = {...
            [-4.5 4.5],...
            [0.25 4.25],...
            [75 325],...
            [10 45],...
            [-4.5 4.5],...
            [0.25 4.25],...
            [75 325],...
            [10 45]}; % part 1
        
        cbarticks = {...
            [-4:2:4],...
            [0 1 2 3 4],...
            [0:50:400],...
            [0:10:40],...
            [-4:2:4],...
            [0 1 2 3 4],...
            [0:50:400],...
            [0:10:40]}; % part 1
        
        cbarminorticks = {...
            [-4:1:4],...
            [0 1 2 3 4],...
            [0:50:400],...
            [0:5:60],...
            [-4:1:4],...
            [0 1 2 3 4],...
            [0:50:400],...
            [0:5:60]}; % part 1
        
        tits = {...
            'T''',...
            '|T''|_{3DST}',...
            '\lambda_{H}',...
            '\lambda_{Z}',...
            'T''',...
            '|T''|_{3DST}',...
            '\lambda_{H}',...
            '\lambda_{Z}'}; % part 1
        
        cbartits = {...
            'K',...
            'K',...
            'km',...
            'km',...
            'K',...
            'K',...
            'km',...
            'km'}; % part 1
        
        cbarlimits = {...
            [-4.5 4.5],...
            [0.25 4.5],...
            [75 325],...
            [5 45],...
            [-4.5 4.5],...
            [0.25 4.5],...
            [75 325],...
            [5 45]}; % usefult to force colorbars to stop
        
        overlay_clims = clims;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Part 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 2
        
        % define isosurfaces:
        iso_cvals = [2 5];
        iso_avals = [0.35 1];
        
        airs_iso_cvals      = [2 5];
        model_iso_cvals     = [1.25 4]; % part 2
        
        
        tbp_isosurfaces = {...
            abs(A.MF2),...
            abs(A.MF),...
            M.MF,...
            M.MF,...
            M.MF,...
            M.MF,...
            M.IN,...
            M.IN}; % part 2
        
        tbp_cvals = {...
            [10 100],...
            [10 100],...
            [10 100],...
            [10 100],...
            airs_iso_cvals,...
            airs_iso_cvals,...
            model_iso_cvals,...
            model_iso_cvals}; % part 2
        
        tbp_avals = {...
            iso_avals,...
            iso_avals,...
            iso_avals,...
            iso_avals,...
            [1 1],...
            [1 1],...
            [1 1],...
            [1 1]}; % part 2
        
        tbp_data = {...
            Airs.MF2,...
            Airs.MF1,...
            Model.MF2,...
            Model.MF1,...
            A.IN,...
            A.IN,...
            M.IN,...
            M.IN}; % part 2
        
        cmaps = {...
            flipud(cbrew('Blues',64)),...
            cbrew('nph_BuOr',64),...
            flipud(cbrew('Blues',64)),...
            cbrew('nph_BuOr',64),...
            cbrew('RdBu',2),...
            cbrew('RdBu',2),...
            cbrew('RdBu',2),...
            cbrew('RdBu',2)}; % part 2
        
        clims = {...
            [-125 0],...
            [-75 25],...
            [-125 0],...
            [-75 25],...
            [-1 1],...
            [-1 1],...
            [-1 1],...
            [-1 1]}; % part 2
        
        cbarticks = {...
            [-500:50:500],...
            [-500:25:500],...
            [-500:50:500],...
            [-500:25:500],...
            [-1 1],...
            [-1 1],...
            [-1 1],...
            [-1 1]}; % part 2
        
        tits = {...
            'MF_X',...
            'MF_Y',...
            'MF_X',...
            'MF_Y',...
            'T''',...
            'T''',...
            'T''',...
            'T'''}; % part 2
        
        cbartits = {...
            'mPa',...
            'mPa',...
            'mPa',...
            'mPa',...
            'K',...
            'K',...
            'K',...
            'K'}; % part 2
        
        cbarlimits = {...
            [-150 0],...
            clims{2},...
            [-150 0],...
            clims{4},...
            [-1 1],...
            [-1 1],...
            [-1 1],...
            [-1 1]}; % usefult to force colorbars to stop at certain values
%         cbarlimits = clims;

        overlay_clims = {...
            [-200 0],...
            clims{2},...
            [-200 0],...
            clims{4},...
            [-1 1],...
            [-1 1],...
            [-1 1],...
            [-1 1]};

end % end switch part.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOTTY PLOTTY PLOT PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch part
    case 1
        axrange = [1 5];
%         axrange = [5 6 7 8];
    case 2
        axrange = [1 2 3 4];
end

for ax = axrange

    axx = subplot(rows,cols,ax);
    
    % set edges to zero:
    tbp_isosurfaces{ax}(edge_indeces(tbp_isosurfaces{ax})) = 0;
    
    % fix an annoying artefact in the model as airs:
    switch part
        case 1
            if ax == 5
                smoo = smoothn(tbp_isosurfaces{ax},[11 11 3]);
                region = {1:20,1:47,1:50};
                tbp_isosurfaces{ax}(region{:}) = smoo(region{:});
            end
            if ax > 4
                smoo = smoothn(tbp_isosurfaces{ax},[15 15 3]);
                region = {1:60,1:9,1:50};
                tbp_isosurfaces{ax}(region{:}) = smoo(region{:});
            end
    end
    
    
    % draw isosurfaces of what ever data on whatever surface:
    for v = 1:length(tbp_cvals{ax})
        
        % PERMUTE EVERYTHING!
        [faces,verteces,colors] = isosurface(D1,D2,D3,permute(tbp_isosurfaces{ax},[2 1 3]),tbp_cvals{ax}(v),permute(tbp_data{ax},[2 1 3]));
%         [faces,verteces,colors] = isosurface(D1,D2,D3,tbp_isosurfaces{ax},tbp_cvals{ax}(v),tbp_data{ax});
        
        patch('Vertices',verteces,'Faces',faces,'FaceVertexCData',colors,...
            'FaceColor','interp','EdgeColor','none','FaceAlpha',tbp_avals{ax}(v),'BackFaceLighting','reverselit');
    
    end
    
    %%% FORMATTING %%%%%%%%%%%%
    
    lighting gouraud
    material dull
    
    % Colors:
    clim([clims{ax}])
    cmap = colormap(gca,cmaps{ax});
    
%     % Map and topography:
%     hMap = surface(xvec,yvec,tpr,Mapr,...
%         'FaceColor','texturemap',...
%         'EdgeColor','none',...
%         'CDataMapping','direct');
    
    % Ticks and Limits:
    axis square
    xlim(xbox); ylim(ybox); zlim([0 80]);
    set(gca,'zgrid','on','ztick',0:20:80)
    set(gca,'xtick',-500:250:500,'ytick',[-400:200:400])
    set(gca,'projection','perspective')
    set(gca,'clipping','on')
    set(gca,'linewi',1.5)
    
    %     %%%%%% Force ZAxis to be front right:
    axx = gca;
    %     axx.ZAxis.FirstCrossoverValue  = axx.XLim(1);
    %     axx.ZAxis.SecondCrossoverValue = axx.XLim(2);
    
    
    % BASIC 3D FORMATTING
    view([-14  20.4])
     
    
    % LIGHTING
    if part == 1
    switch ax
        case {1,5}
            camlight(-90,45); % creates nice shadows
            c1 = camlight('left');
        otherwise
            c2 = camlight('right');
            c3 = camlight('headlight');
    end
    else
        % For part 2, just do all of them from the left
        camlight(-90,45); % creates nice shadows
        c1 = camlight('left');
        c3 = camlight('headlight');
    end
    
    % AXES TITLES
    if ax < 5
        nph_text([0.5 1.05],tits{ax},'fontsize',1.5*fs,'fontweight','bold');
    end
    
    % Pivot Colorbars for Meridional MF:
    if part == 2
        switch ax
            case {2,4}
                pivval = 0;
                cl = clims{ax};
                cl = cl - clims{ax}(1);
                pivfrac = abs((pivval-clims{ax}(1)) ./ cl(2));
                cmaps{ax} = cbrew('nph_BuOr',64,'pivot',pivfrac);
                colormap(gca,cmaps{ax});
            otherwise
                pivfrac = 0.5;
        end
    end
    
    % COLORBARS
    if part == 1
        cbar_axrange = {5,6,7,8};
    else
        cbar_axrange = {1,2,3,4};
    end
    
    switch ax
        case cbar_axrange
            
            % COLORBARS
            cbar = colorbar('location','southoutside');
            
            apos = get(gca,'Position');
            cbar.Position = apos;
            
            % reposition
            wid = 0.55;
            cbar.Position(1) = apos(1) + ((1-wid)/2)*apos(3);
            cbar.Position(2) = apos(2) - 0.2*apos(4);
            cbar.Position(3) = wid   * apos(3);
            cbar.Position(4) = 0.06 * apos(4);
            
            cbar.Ticks          = cbarticks{ax};
            cbar.Label.String   = cbartits{ax};
            cbar.Limits         = clims{ax};
            cbar.FontSize       = fs;
            cbar.LineWidth      = 1.5;
            cbar.TickDirection  = 'out';
            
            % Force limits:
            cbar.Limits = cbarlimits{ax};
            
            % minor ticks?
            switch ax
                case {5,8}
                cbar.Ruler.MinorTick = 'on';
                cbar.Ruler.MinorTickValues = cbarminorticks{ax};
            end
            
    end
    
    % X and y labels - xlabel etc are shit for 3D
    if part == 1
        xlabel_axrange = {1,2,3,4,5,6,7,8};
        ylabel_axrange = {1,5};
    else % part 2
        xlabel_axrange = {1,2,3,4};
        ylabel_axrange = {1};
    end
    switch ax
        case xlabel_axrange
            xl = xlabel('_{x (km)}');
            xl.Position = [0 -750 0];
    end
    switch ax
        case ylabel_axrange
            yl = ylabel('_{y (km)}');
            yl.Position = [-1000 0 0];
    end
    
    % Z Labels
    if ax == 1 || ax == 5
        zl = zlabel('Altitude (km)');
    end
    if ax == 5
        zl = zlabel('Altitude (km)');
    end
    
    
    % Any more axes specific formatting....?
    if part == 1
    switch ax
        case 1
            % AIRS/MODEL LABELS
            text(-1400,0,55,'AIRS','fontsize',2*fs,'rotation',90,'horizontalalignment','center');
        case 2
        case 3
        case 4
        case 5
            % AIRS/MODEL LABELS
            text(-1400,0,55,'Model as AIRS','fontsize',2*fs,'rotation',90,'horizontalalignment','center');
        case 6
    end
    end
    
    %%%%%% MODEL AND AIRS DOMAINS --------------------------------------
    ModelXLim = [-700 700];
    ModelYLim = [-550 550];
    if part == 1
        airs_domain_axrange = {1,2,3,4};
    else
        airs_domain_axrange = {1,2};
    end
    switch ax
        case airs_domain_axrange
            % AIRS DOMAIN
            axx = gca;
            % top
            hold on; plot3(ModelXLim([1 2 2 1 1]),ModelYLim([1 1 2 2 1]),60*ones(1,5),'color',[1 1 1 0.75],'linewi',1.5)
            hold on; plot3(ModelXLim([1 2 2 1 1]),ModelYLim([1 1 2 2 1]),60*ones(1,5),'color',[c_color('red') 0.5],'linewi',1.5,'linest','--')
            % bottom
            hold on; plot3(ModelXLim([1 2 2 1 1]),ModelYLim([1 1 2 2 1]),20*ones(1,5),'color',[1 1 1 0.75],'linewi',1.5)
            hold on; plot3(ModelXLim([1 2 2 1 1]),ModelYLim([1 1 2 2 1]),20*ones(1,5),'color',[c_color('red') 0.5],'linewi',1.5,'linest','--')
        otherwise
            % MODEL DOMAIN
            axx = gca;
            % top
            hold on; plot3(ModelXLim([1 2 2 1 1]),ModelYLim([1 1 2 2 1]),75*ones(1,5),'color',[1 1 1 0.75],'linewi',1.5)
            hold on; plot3(ModelXLim([1 2 2 1 1]),ModelYLim([1 1 2 2 1]),75*ones(1,5),'color',[c_color('blue') 0.5],'linewi',1.5,'linest','--')
            % bottom
            hold on; plot3(ModelXLim([1 2 2 1 1]),ModelYLim([1 1 2 2 1]),1.5*ones(1,5),'color',[1 1 1 0.75],'linewi',1.5)
            hold on; plot3(ModelXLim([1 2 2 1 1]),ModelYLim([1 1 2 2 1]),1.5*ones(1,5),'color',[c_color('blue') 0.5],'linewi',1.5,'linest','--')
    end

    
   % AXES LABEL
    nph_text([0.1 0.95],['(' alphabet(ax) ')'],'fontsize',1.75*fs,'horizontalalignment','left');
    
    %%% NORTH ARROW
    arrowlength = 200;
    arrowstart = [-550 -425];
    arrowcolor = [.95 .95 .95];
    hold on; text(arrowstart(1),arrowstart(2)+arrowlength+50,1,'N','color',arrowcolor,'rotation',20,'fontsize',0.9*fs,'fontweight','bold','HorizontalAlignment','center','VerticalAlignment','middle')
    hold on; quiver3(arrowstart(1),arrowstart(2),1,0,arrowlength,0,'color',arrowcolor,'linewi',1.5,'maxheadsize',0.8)
    cross_x = arrowstart(1) + [-20 20];
    cross_y = arrowstart(2) + 0.4*arrowlength .* [1 1];
    hold on; plot3(cross_x,cross_y,[1 1],'color',arrowcolor,'linewi',1.5)
    
    
    % FONTS
    setfont(fs)
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% OVERLAY SLICE AT 40KM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% overlay a slice at 40km on top of each plot. Set regions outside the
% biggest isoblob to be invisible.

if overlay

    zlev = 27; % 40.5km
    
    axx = gca;
    axx2 = axes('position',axx.Position,'xlim',axx.XLim,'ylim',axx.YLim,'zlim',axx.ZLim,...
        'color','none','xcolor','none','ycolor','none','zcolor','none',...
        'projection','perspective','clipping','on','view',axx.View);
    
    axis square
    
    xmid = -275; xwid = 400;
    ymid = -100; ywid = 325;
    [X,Y] = meshgrid(linspace(xmid-xwid,xmid+xwid,120),linspace(ymid-ywid,ymid+ywid,90));
    Z = 80*ones(size(X));
    
    overlay_data = tbp_data{ax}(:,:,zlev);
    
    if part == 1
        if ax == 1 || ax == 5
            overlay_cmap = cbrew('RdBu',31);
        else
            overlay_cmap = cmaps{ax};
            overlay_isodata = tbp_isosurfaces{ax}(:,:,zlev);
            overlay_data(smoothn(overlay_isodata,[5 5]) < 0.8) = NaN;
        end
    else % part 2
            overlay_cmap = cmaps{ax};
            overlay_isodata = tbp_isosurfaces{ax}(:,:,zlev);
            overlay_data = smoothn(tbp_data{ax}(:,:,zlev),[5 5]);
            overlay_data(smoothn(overlay_isodata,[5 5]) < 5) = NaN;
    end
    
    
    hold on; S = surf(X,Y,Z,overlay_data);
    S.LineStyle = 'none';
    
    colormap(axx2,overlay_cmap);
    clim(overlay_clims{ax})
    
    % line around the edge:
    hold on; plot3(X(edge_indeces(X)),Y(edge_indeces(Y)),Z(edge_indeces(Z)),'color','k','linewi',1.5);
    
    % coastline?
    LonBox = [-50 -30];
    LatBox = [-60 -50];
    C = nph_draw_coastline([LonBox(1) LatBox(1) ; LonBox(2) LatBox(2)],0,0,'noplot','color','k');
    for i = 1:length(C)
        if length(C(i).Lon) > 1
            [d,az] = distance(clat,clon,C(i).Lat,C(i).Lon);
            xx = deg2km(d).*sind(az);
            yy = deg2km(d).*cosd(az);
            zz = 80 * ones(size(d));
            % modify for shrinking:
            xx = xx .* 0.6;
            yy = yy .* 0.6;
            hold on; plot3(xx+xmid,yy+ymid,zz,'color','w','linewi',1.5);
            hold on; plot3(xx+xmid,yy+ymid,zz,'color','k','linewi',1);
        end
    end
    
    
end % end if overlay    


    



    drawnow;
    
    
    
end % next axes





return




%% EXPORT? ================================================================

switch part
    case 1
        savename = ['~/Desktop/' granule '_AirsModel3D_MultiProp_pt1'];
    case 2
        savename = ['~/Desktop/' granule '_AirsModel3D_MultiProp_pt2'];
end

disp(['Exporting to "' savename '"...'])
nph_saveas(gcf,savename,'png')

disp('Done.')



return










































