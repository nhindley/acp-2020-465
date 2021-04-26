

% Plot the baseline geography for the SGWEX Model domain. Or any other
% domain for that matter, though you'll need to switch the imagery map to a
% global one (i made a SAAP region one to speed up loading)

function hMap = sgwex_plot_baseline_geography(LonBox,LatBox,varargin)

matlabdirec = '/Users/neil/Drive/MATLAB/';

if nargin == 0
    % outer box region:
    LonBox = [-47.5 -26.5];
    LatBox = [-59 -49];
end

%define region to work in
minlon = LonBox(1); maxlon = LonBox(2);
minlat = LatBox(1); maxlat = LatBox(2);

%% LOAD TOPOGRAPHY AND MAP IF REQUIRED ====================================
% See if Map and Topography already loaded in base workspace:
try
    tp = evalin('base','tp');
    Map = evalin('base','Map');
catch % otherwise load:
    disp('Topography loading...')
    tp = load([matlabdirec 'topography/easy_tenth_degree_topography/Global.mat']);
    % gridded interpolant:
    Ftp = griddedInterpolant({-90:0.1:90,-180:0.1:180},tp.tp,'linear','none');
    % regionally subset
    rlons = minlon:0.05:maxlon; % over interpolate for smoothness
    rlats = minlat:0.05:maxlat;
    % evaulate
    tp.elev = Ftp({rlats,rlons});
    [tp.lons,tp.lats] = meshgrid(rlons,rlats);
    tp.elev(tp.elev < 0) = 0; %shallower sea
    tp.elev = tp.elev ./ 1000; % CONVERT TO KM!!!
    
    % %% LOAD MAP
    disp('Map Loading...') % each degree is 60 pixels, i think!
    
    load([matlabdirec 'imagery/SAAP_Region_Map.mat'])
    FMap = griddedInterpolant({Map.Lat,Map.Lon,1:3},single(Map.Map),'linear','none');
    % regionally subset
    rlons = minlon:0.01:maxlon; % over interpolate for smoothness
    rlats = minlat:0.01:maxlat;
    Map.Map = uint8(FMap({rlats,rlons,1:3}));
end

%% PLOT BASE GEOGRAPHY ====================================================

%%% MAP ---------------------------------------------------------------
tp.elev(tp.elev == 0) = 0;
hMap = surface(tp.lons,tp.lats,2*tp.elev,Map.Map,...
    'FaceColor','texturemap',...
    'EdgeColor','none',...
    'CDataMapping','direct');

%%% COASTLINE ---------------------------------------------------------
coastlinespec = {'color',c_color('yellow'),'linewi',1};
nph_draw_coastline([LonBox(1) LatBox(1) ; LonBox(2) LatBox(2)],0.25,coastlinespec{:});

%%% POSITIONING -------------------------------------------------------

xlim(LonBox); % Lon
ylim(LatBox); % Lat
zlim([0 80]); % Alt

camview = [-4.265 18];

set(gca,'view',camview);

%%% SUBPLOT CONDITIONAL FORMATTING ------------------------------------

set(gca,'zgrid','on');

set(gca,'linewi',2,'xtick',LonBox,'ytick',LatBox);
set(gca,'xgrid','off','ygrid','off');
set(gca,'clipping','off');

lonrange = -45:5:-30;
latrange = -55:5:-50;

set(gca,'xtick',lonrange,'xticklabels',lonlabels(lonrange));
set(gca,'ytick',latrange,'yticklabels',latlabels(latrange));
set(gca,'ztick',0:20:80);


setfont(30);

%%
%%% OTHER AXES STUFF:
if ~isempty(varargin)
    set(gca,varargin{:});
end




end
