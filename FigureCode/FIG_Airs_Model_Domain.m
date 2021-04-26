

%%

% Plot opening figure of Airs vs SGWEX Model paper, a map of the island with
% the AIRS and Model domains:


%% DIRECTORIES ============================================================

% modeltype = 'Model';
modeltype = 'Model_as_AIRS';

matlabdirec = '/Users/neil/Drive/MATLAB/';
airsdirec = '/Users/neil/data/sgwex/3DST/Airs/';
% modeldirec = '/Users/neil/data/sgwex/3DST/Model/downsized/';
modeldirec = '/Users/neil/data/sgwex/3DST/Model_as_AIRS/downsized/';




% %% CREATE GRID ============================================================
%
% sgbox = sgwexbox;
%
% % at set radiius fro SG
%
% dr = 1; % km
% r = 0:dr:1000; % km
%
% dt = 0.5;
% theta = 0:dt:360;
%
% [R,Th] =  meshgrid(r,theta);
%
% [LAT,LON] = reckon(sgbox.centrelat,sgbox.centrelon,km2deg(R),Th);
%
% %%
%
% load([matlabdirec 'imagery/map.mat'])
% % load([matlabdirec 'imagery/SAAP_Region_Map.mat'])
%
% % distance from lower left corner?
% [MLat,MLon] = meshgrid(Map.Lon,Map.Lat);
%
%
% [arclen,az] = distance(sgbox.centrelat,sgbox.centrelon,MLat,MLon);
%     xdist = deg2km(arclen) .* sind(az);
%     ydist = deg2km(arclen) .* cosd(az);
%
% FMap = griddedInterpolant({xdist,ydist,1:3},single(Map.Map),'linear','none');
% % FMap = griddedInterpolant({Map.Lat,Map.Lon,1:3},single(Map.Map),'linear','none');
%
% % Map.Map = uint8(FMap(repmat(LAT,1,1,3),repmat(LON,1,1,3),permute(repmat(1:3))));







%% LOAD MAP AND TOPOGRAPHY ================================================
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sgbox = sgwexbox;

% outer box region:
% LonBox = [-47.5 -26.5];
% LatBox = [-59 -49];

LonBox = sgbox.centrelon + plusminus(40);
LatBox = sgbox.centrelat + [-30 20];

%define region to work in
minlon = LonBox(1); maxlon = LonBox(2);
minlat = LatBox(1); maxlat = LatBox(2);

% %% LOAD TOPOGRAPHY ======================================================
disp('Topography loading...')
tp = load([matlabdirec 'topography/easy_tenth_degree_topography/Global.mat']);

% gridded interpolant:
Ftp = griddedInterpolant({-90:0.1:90,-180:0.1:180},tp.tp,'linear','none');

% regionally subset
rlons = minlon:0.05:maxlon; % over interpolate for smoothness
rlats = minlat:0.05:maxlat;

% evaulate:
tp.elev = Ftp({rlats,rlons});
[tp.lons,tp.lats] = meshgrid(rlons,rlats);

tp.elev(tp.elev < 0) = 0; %shallower sea
tp.elev = tp.elev ./ 1000; % CONVERT TO KM!!!

% %% LOAD MAP =============================================================
disp('Map Loading...')
% each degree is 60 pixels, i think!

% load([matlabdirec 'imagery/SAAP_Region_Map.mat'])
load([matlabdirec 'imagery/map.mat'])
FMap = griddedInterpolant({Map.Lat,Map.Lon,1:3},single(Map.Map),'linear','none');

% regionally subset
rlons = minlon:0.01:maxlon; % over interpolate for smoothness
rlats = minlat:0.01:maxlat;

Map.Map = uint8(FMap({rlats,rlons,1:3}));

% Brighten map:
Map.Map = Map.Map * 1.2;
Map.Map(Map.Map > 255) = 255;


%% "DISTANCE FROM SOUTH GEORGIA" GRID =======================================

[arclen,az] = distance(sgbox.centrelat,sgbox.centrelon,tp.lats,tp.lons);
xdist = deg2km(arclen) .* sind(az);
ydist = deg2km(arclen) .* cosd(az);

Re = 6371;
zdist = Re - (Re .* sind(arclen));
    
figure; hold all; grid on; axis square;

hMap = surface(xdist,ydist,2*tp.elev,Map.Map,...
    'FaceColor','texturemap',...
    'EdgeColor','none',...
    'CDataMapping','direct');

set(gca,'clipping','off');
set(gca,'cameratarget',[0 0 0]);

% zlim([0 100])



return









%% PLOT BASE GEOGRAPHY ====================================================

figure; hold all; grid on;
whitefig;

% figpos([168         212        1631         753])
%
% vert_gap = 0;       horz_gap = 0.1;
% lower_marg = 0.1;   upper_marg = 0.125;
% left_marg = 0.08;   right_marg = 0.08;
%
% rows = 1; cols = 2;
%
% subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

sgwex_plot_baseline_geography(LonBox,LatBox);
set(gca,'clipping','on');
drawnow;


zlabel('Altitude (km)')
sgwex_plot_airs_model_box('model');








return




[sgx,sgy,sgz] = lla2ecef(sgbox.centrelat,sgbox.centrelon,(6371*1000));


% [X,Y,Z] = lla2ecef(PLAT,PLON,(tp.elev*1000)+(6371*1000));
[X,Y,Z] = lla2ecef(tp.lats,tp.lons,(tp.elev*1000)+(6371*1000));

figure; hold all; grid on; axis square;

hMap = surface(X,Y,Z,Map.Map,...
    'FaceColor','texturemap',...
    'EdgeColor','none',...
    'CDataMapping','direct');

set(gca,'cameratarget',[sgx sgy sgz]);


% zlim([0 6371*1000])


return






%%

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


whitefig; setfont(30);




return



%% PLOT BASE GEOGRAPHY ====================================================

figure; hold all; grid on;
whitefig;

figpos([168         212        1631         753])

vert_gap = 0;       horz_gap = 0.1;
lower_marg = 0.1;   upper_marg = 0.125;
left_marg = 0.08;   right_marg = 0.08;

rows = 1; cols = 2;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

% plot baseline geography, Boxes and titles ===============================
for ax = 1:2
    subplot(rows,cols,ax);
    sgwex_plot_baseline_geography(LonBox,LatBox); drawnow;
    switch ax
        case 1 % Airs
            ax1 = gca;
            %             axis vis3d
            zlabel('Altitude (km)')
            sgwex_plot_airs_model_box('airs');
            title({'AIRS',' '})
            
        case 2 % Model as AIRS
            ax2 = gca;
            %             axis vis3d
            sgwex_plot_airs_model_box('model');
            title({'Model at AIRS resolution',' '})
            
    end
end



%%















%% Find 3D AIRS granules (or sequential granule pairs) that completely
% overlap the SGWEX modelling domain box.



% % % % % % % SGWEX BOX:
% % % % % minlon = -42.5; maxlon = -31.71;
% % % % % minlat = -58.55; maxlat = -50.46;
% % % %
% % % % % WARPING-CORRECTED SGWEX BOX:
% % % % % True centre location of box.
% % % % centrelon = -37.1067443847656;
% % % % centrelat = -54.5067496299744;
% % % %
% % % % % Box width ~1200x900, diagonals = 1500km
% % % % % 4 --- 3
% % % % % |     |
% % % % % 1 --- 2
% % % %
% % % % % NOTE: don't use the referenceEllipsoid/referenceSphere functions, I can't
% % % % % seem to get a correct answer with them and i don't know why. Had this
% % % % % trouble at Dstl too.
% % % %
% % % % [boxlat(1),boxlon(1)] = reckon(centrelat,centrelon,km2deg(750),180+45);
% % % % [boxlat(2),boxlon(2)] = reckon(centrelat,centrelon,km2deg(750),180-45);
% % % % [boxlat(3),boxlon(3)] = reckon(centrelat,centrelon,km2deg(750),+45);
% % % % [boxlat(4),boxlon(4)] = reckon(centrelat,centrelon,km2deg(750),-45);
% % % %
% % % %
% % % % % boxlats = [boxlat boxlat(1)];
% % % % % boxlons = [boxlon boxlon(1)];

figure; hold all; grid on; whitefig;
figpos([591 165 1144 880])


sgwexbox;

m_proj('ortho','lat',centrelat,'lon',centrelon,'radius',35)
m_coast('patch',c_color('light_grey'));
m_grid('xtick',[],'ytick',[],'color',get(gcf,'color'));
setfont(20)


% Allow for granules that just miss out longitudinally:
tolerance = 1;
boxlon([1 4]) = boxlon([1 4]) + tolerance;
boxlon([2 3]) = boxlon([2 3]) - tolerance;



BoundingBox = [-80 -80 ; 0 0];

nph_draw_coastline(BoundingBox,0,'color','k','linewi',1,'mplot');

% hold on; m_plot(boxlon([1:end 1]),boxlat([1:end 1]),'color','r','linewi',2);

% On selecting granules: Having been looking around the model box and AIRS
% swaths for some time now, there are generally two times when AIRS does
% an overpass that looks good over South Georgia, which sits on the
% interface between two granules divisions, generally. These occur around
% granules 031-032 and 166-167, or thereabouts due to precession of the
% orbits. It's quite rare that a single granule covers all the box.
%
% Typically, an overpass might have three corners of the modelling box
% covered, at least, with the island itself always being contained within
% the granule. So I think it's safe to take these granules that miss
% sections of the box and 3DST them anyway with zeros for missing values.
%

%% LOAD AIRS 3DST =========================================================

direc = '/Volumes/:/Users/neil/data/AIRS/airs3d/SGWEX-SGeorgia/';
files = dir([direc 'airs*.nc']);

savedirec = '/Volumes/:/Users/neil/data/AIRS/airs3d/SGWEX-SGeorgia/three_corner_overlapping_granule_pairs/';

% sort files into ascending order:


for i = 1:length(files)-1
    % for i = 1:100,
    
    % Select sequential pairs of granules:
    str1 = strsplit(files(i).name,{'_','.'});
    currentnum = str2double(str1{4});
    
    str2 = strsplit(files(i+1).name,{'_','.'});
    nextnum = str2double(str2{4});
    
    if nextnum ~= currentnum+1
        %         disp('not a pair')
        continue
        %     else
        %         disp({str1{3},str1{4},str2{4}})
    end
    
    % Criteria for granule pair selection:
    % - separately, both granules must each overlap at least one corner of the box.
    % - together, both granules must overlap at least three corners.
    % This places the granule interface inside the sgwex box, and kind of
    % centres the granules over the island.
    
    nc1 = getnet([direc files(i).name]);
    nc2 = getnet([direc files(i+1).name]);
    %%
    nc1.ret_temp = permute(nc1.ret_temp,[2 3 1]);
    nc2.ret_temp = permute(nc2.ret_temp,[2 3 1]);
    
    Airs = struct;
    
    Airs.T      = cat(2,nc1.ret_temp,nc2.ret_temp);
    Airs.Lon    = cat(2,  nc1.l1_lon,  nc2.l1_lon);
    Airs.Lat    = cat(2,  nc1.l1_lat,  nc2.l1_lat);
    Airs.Alt = nc2.ret_z;
    
    % Lars' Airs time is seconds since 2000?
    starttimeindayssince2000 = nc1.l1_time(1,1)/60/60/24;
    Airs.StartTime = datenum('01-Jan-2000') + starttimeindayssince2000;
    
    % granule 1 coords
    lons1 = [nc1.l1_lon(1:90,1)' nc1.l1_lon(90,1:135) ...
        nc1.l1_lon(90:-1:1,135)' nc1.l1_lon(1,135:-1:1)];
    lats1 = [nc1.l1_lat(1:90,1)' nc1.l1_lat(90,1:135) ...
        nc1.l1_lat(90:-1:1,135)' nc1.l1_lat(1,135:-1:1)];
    
    % granule 2 coords
    lons2 = [nc2.l1_lon(1:90,1)' nc2.l1_lon(90,1:135) ...
        nc2.l1_lon(90:-1:1,135)' nc2.l1_lon(1,135:-1:1)];
    lats2 = [nc2.l1_lat(1:90,1)' nc2.l1_lat(90,1:135) ...
        nc2.l1_lat(90:-1:1,135)' nc2.l1_lat(1,135:-1:1)];
    
    % granule pair coords
    pairlons = [Airs.Lon(1:end,1)' Airs.Lon(end,1:end) ...
        Airs.Lon(end:-1:1,end)' Airs.Lon(1,end:-1:1)];
    pairlats = [Airs.Lat(1:end,1)' Airs.Lat(end,1:end) ...
        Airs.Lat(end:-1:1,end)' Airs.Lat(1,end:-1:1)];
    %%
    if any(inpolygon(boxlon,boxlat,lons1,lats1)) ...
            && any(inpolygon(boxlon,boxlat,lons2,lats2))
        
        if sum(inpolygon(boxlon,boxlat,pairlons,pairlats)) >= 3
            
            % Plot? with interface line shown in dots
            hold on; m_plot(lons1,lats1,'color','k','linestyle',':','linewi',0.5);
            hold on; m_plot(pairlons,pairlats,'color','k');
            
            drawnow;
            
            % Handle time for the save string:
            % We bin overpasses into the 0300 and 1700 overpasses
            
            hours = round((Airs.StartTime - floor(Airs.StartTime)) * 24);
            
            if hours >= 2 && hours <= 4
                hours = 3;
            end
            if hours >= 16 && hours <= 18
                hours = 17;
            end
            
            % Use datestr to handle zeros and shit:
            %             roundedtime = floor(Airs.StartTime) + (hours/24);
            %
            %             str = datestr(roundedtime,'yyyymmdd-HHMM');
            %
            %             savefilename = [str '_AirsSG.mat'];
            %
            %             save([savedirec savefilename],'Airs');
            %
            %             disp(['Saved "' savefilename '"'])
            
        end
        
    end
    
end

% %%
% MODEL BOX:
hold on; m_plot(Model.Lon(edge_indeces(Model.Lon)),Model.Lat(edge_indeces(Model.Lat)),'color','r','linewi',2);



return


%% EXPORT

nph_saveas(gcf,'~/Desktop/overlapping_granules','png')





























% Find out if they overlap the box:
%     (at least three corners overlapping)

%     if sum(inpolygon(boxlon,boxlat,airslons,airslats)) >= 3
%
%         hold on; plot(airslons,airslats,'color','k');
%         drawnow;
%
% %     else
%
%         continue
%
%         title(strrep(files(i).name,'_','-'));
%
%         drawnow;
%
%         disp(files(i).name)
%
%         % save([direc '../overlapping_swaths/fullyoverlapping/' files(i).name],'Airs')
%
%         savestr = [direc '../overlapping_swaths/fullyoverlapping/' currentstr(1:end-4) '-' nextstr(end-6:end-4) '.mat'];
%
%         save(savestr,'Airs');
%
%         disp(['Saved granules ' currentstr ' and ' nextstr])
%
%
%
% %     end
%
%
% end

%         % If it's a sequential pair, do stuff:
%         if nextnum == currentnum + 1,
%
%             T1 = Airs.T;
%             Lat1 = Airs.Lat;
%             Lon1 = Airs.Lon;
%             alt = Airs.Alt;
%
%             load([direc nextstr])
%
%             T2 = Airs.T;
%             Lat2 = Airs.Lat;
%             Lon2 = Airs.Lon;
%
%             % Concatenate in along-track direction:
%
%             Airs = struct;
%
%             Airs.T = cat(2,T1,T2);
%             Airs.Lat = cat(2,Lat1,Lat2);
%             Airs.Lon = cat(2,Lon1,Lon2);
%             Airs.Alt = alt;
%
%             %             T1 = Airs.T;
%             %             Lat1 = Airs.Lat;
%             %             Lon1 = Airs.Lon;
%             %
%             %             Airs = struct;
%             %             Airs.T = T1;
%             %             Airs.Lat = Lat1;
%             %             Airs.Lon = Lon1;
%
%
%             sz = size(Airs.T);
%
%             x = [1:sz(1) sz(1)*ones(1,sz(2)) sz(1):-1:1 ones(1,sz(2)) 1];
%             y = [ones(1,sz(1)) 1:sz(2) sz(2)*ones(1,sz(1)) sz(2):-1:1 1];
%
%             airslons = [Airs.Lon(1:sz(1),1)' Airs.Lon(sz(1),1:sz(2)) ...
%                 Airs.Lon(sz(1):-1:1,sz(2))' Airs.Lon(1,sz(2):-1:1)];
%             airslats = [Airs.Lat(1:sz(1),1)' Airs.Lat(sz(1),1:sz(2)) ...
%                 Airs.Lat(sz(1):-1:1,sz(2))' Airs.Lat(1,sz(2):-1:1)];

%     if all(inpolygon([minlon minlon maxlon maxlon],[minlat maxlat maxlat minlat],airslons,airslats)),
%             if all(inpolygon(boxlon,boxlat,airslons,airslats)),
%
%                 hold on; plot(airslons,airslats,'color','k');
%
%                 title(strrep(files(i).name,'_','-'));
%
%                 drawnow;
%
%                 disp(files(i).name)
%
%                 % save([direc '../overlapping_swaths/fullyoverlapping/' files(i).name],'Airs')
%
%                 savestr = [direc '../overlapping_swaths/fullyoverlapping/' currentstr(1:end-4) '-' nextstr(end-6:end-4) '.mat'];
%
%                 save(savestr,'Airs');
%
%                 disp(['Saved granules ' currentstr ' and ' nextstr])
%
%
%
%             end
%
% end
%     end

%     x = [1:sz(1) sz(1)*ones(1,sz(2)) sz(1):-1:1 ones(1,sz(2)) 1];
%     y = [ones(1,sz(1)) 1:sz(2) sz(2)*ones(1,sz(1)) sz(2):-1:1 1];
%
%     lons = [Airs.Lon(1:sz(1),1)' Airs.Lon(sz(1),1:sz(2)) ...
%         Airs.Lon(sz(1):-1:1,sz(2))' Airs.Lon(1,sz(2):-1:1)];
%     lats = [Airs.Lat(1:sz(1),1)' Airs.Lat(sz(1),1:sz(2)) ...
%         Airs.Lat(sz(1):-1:1,sz(2))' Airs.Lat(1,sz(2):-1:1)];
%
%
%     % Check if within SGWEX box:
%
%     if all(inpolygon([minlon minlon maxlon maxlon],[minlat maxlat maxlat minlat],lons,lats)),
%
%         hold on; plot(lons,lats,'color','k');
%
%         title(strrep(files(i).name,'_','-'));
%
%         disp(files(i).name);
%
%         drawnow;
%
%     end
%
%
%
% end
%






















