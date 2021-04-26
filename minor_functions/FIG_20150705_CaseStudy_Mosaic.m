

% In response to reviewer #3, they didn't believe that you could get GWs at
% sort horizontal scales with massive perturbations. They thought that
% these weren't real and instead were due to the "coarse" vertical grid in
% the model.

% Well, on 20150705-1700, the AIRS measurements are aligned perfectly for a
% case study. The madir fall directly over the island, and the across-track
% direction is aligned exactly opposite to the orinetation of the mountain
% wave field propagating into the mean wind.

% this provides the best possible AIRS horizontal resolution over the island.

% by re-gridding the model on the AIRS grid, we can see how similar the
% AIRS and model really were for this example.

% also take a cut along the across-track direction to show the size and
% horizontal scale of these perturbations.

% 2x2 panel of M.T, M.Tp, A.T, A.Tp, where the latter are all on the airs
% grid this time, not the model grid.

% underneath, a line plot of a cross-track row over the island.




%%%% DIRECTORIES ============================================================

airsdirec           = '/Volumes/SDBlue/data/sgwex/AIRS/three_corner_overlapping_granule_pairs/matlab/';
modeldirec          = '/Volumes/SDBlue/data/sgwex/3DST/Model/';

%%%% LOAD TOPO AND MAP ======================================================
% need this for the amp/gini/distance plots :)

disp('Loading Topography...')
load('/Users/neil/Drive/MATLAB/topography/SRTM/srtm_29_23/srtm_29_23.mat')

% trim to island area
topo = double(Topo.Topo(1:1500,2000:5500));
topolat = Topo.Lat(1:1500);
topolon = Topo.Lon(2000:5500);

% fix sea level
topo(topo < 0) = 0;

% % smooth with a gaussian to simulate AIRS footprint
% spacing = [mean(diff(topolat)) mean(diff(topolon))];
% fwhm = [0.0125 0.03]; % degrees, guessing...
% topo_sm = imgaussfilt(topo,(fwhm./2.355)./spacing,'padding','replicate');
%
% disp('Creating topography interpolant...')
%
% 
% 
% % make an interpolant for later:
% [x,y] = meshgrid(topolon,topolat);
% F_topo = scatteredInterpolant(x(:),y(:),topo_sm(:),'linear','none');
% 


%
% clat = -54.5269874182436;
% clon = -37.1367495437688;
%
% [d,az] = distance(clat,clon,repmat(clat,1,length(Topo.Lon)),Topo.Lon);
% xdist = deg2km(d).*sind(az);

% % smooooth
% topo = movmean(topo,51);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LOAD AIRS AND MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

granule = '20150705-1700';

disp('Loading AIRS and full-resolution model...')

% AIRS
load([airsdirec granule '_AirsSG.mat'])

% full resolution model wind:
load('/Users/neil/data/20150705-1700_SG_uvwtpd.mat')
modelwind.u = permute(Model.u,[2 1 3]);
modelwind.v = permute(Model.v(:,1:600,:),[2 1 3]);
modelwind.w = permute(Model.w,[2 1 3]);
modelwind.z = Model.Alt_1;

% Full-resolution Model
load([modeldirec granule '_Model_3DST_4dp.mat'])





% =========================================================================
% EXTRACT AIRS LAYER

zlev = 16; % 45km

A.Tg = Airs.T;
A.Tg_lev = A.Tg(:,:,zlev);

[A.Tpg_lev,A.bg_lev] = nph_airs_4dp_detrend(A.Tg_lev,1,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PUT MODEL ON AIRS GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Putting model on AIRS grid...')

M.T = double(Model.T);

% extract level
[~,zind] = min((Model.Alt - Airs.Alt(zlev)).^2);
M.T_lev = M.T(:,:,zind);

% % % % % MODEL ON AIRS FOOTPRINTS - not used now
% % % % % grid_spacing = [1.5 1.5];
% % % % % sig = ([13.5 13.5] ./ grid_spacing) ./ 2.355; % see hoffmann 2014, footprint size = 13.5km
% % % % % M.T_sm = imgaussfilt(M.T,sig,'padding','replicate');

% MODEL ON AIRS FOOTRPINTS AND AT AIRS VERT RES
% vert_res = airs3d_vert_res(Airs.Alt(zlev));
vert_res = 10; % 10 km, for z = 45km

grid_spacing = [1.5 1.5 1.5];
sig = ([13.5 13.5 vert_res] ./ grid_spacing) ./ 2.355; % see hoffmann 2014, footprint size = 13.5km
MA.T_sm = imgaussfilt3(M.T,sig,'padding','replicate');

% put on airs grid:
F_ma = scatteredInterpolant(Model.Lon(:),Model.Lat(:),linearise(MA.T_sm(:,:,zind)),'linear','none');
MA.Tg_lev = F_ma(Airs.Lon,Airs.Lat);

% add airs retrieval noise:
disp('adding AIRS noise to the model...')
% 20150620-1700 looks good.
% best to use a real granule cos the noise is different for each altitude.
B = load(['/Volumes/SDBlue/data/sgwex/AIRS/three_corner_overlapping_granule_pairs/matlab/20150620-1700_AirsSG.mat']);
[B.Airs.tp,~] = nph_airs_4dp_detrend(B.Airs.T,1,4);
airsnoise = B.Airs.tp(:,136:270,:);
airsnoise = cat(2,airsnoise,airsnoise);
% randomise each layer to get rid of the join:
asz = size(B.Airs.T);
for z = 1:asz(3)
    layer = airsnoise(:,:,z);
    layer = reshape(layer(randperm(prod(asz(1:2)))),asz(1:2));
    airsnoise(:,:,z) = layer;
end
% % put the airs noise on the regular vertical grid:
% FF = griddedInterpolant({1:90,1:270,A.Airs.Alt},airsnoise,'linear','none');
% airsnoise = FF({1:90,1:270,Model.Alt});
% and add it to the model as airs:
MA.Tg_lev = MA.Tg_lev + (1*airsnoise(:,:,zlev));
% and you're done!

% full resolution model interpolant for XSECS later
F_m = scatteredInterpolant(Model.Lon(:),Model.Lat(:),linearise(M.T(:,:,zind)),'linear','none');

% get perts on cross track grid
[MA.Tpg_lev,MA.bgg_lev] = nph_airs_4dp_detrend(MA.Tg_lev,1,4);

% also get full resolution model wind for "inspection"...
% need a new levelf for this:
[~,windind] = min(((modelwind.z./1000) - Airs.Alt(zlev)).^2);
F_wind = scatteredInterpolant(Model.Lon(:),Model.Lat(:),linearise(double(modelwind.u(:,:,windind))),'linear','none');


% 
% % also, for the vertical profile, put the model onto the AIRS grid and
% % apply the vertical resolution for each height
% MA.T_sm = nan(size(M.T));
% vr = airs3d_vert_res(Model.Alt);
% uvr = unique(vr);
% for v = 1:length(uvr)
%     disp(uvr(v))
%     % locations for this vert res:
%     zinds = uvr(v) == vr;
%     % smooth away
%     sig = ([13.5 13.5 uvr(v)] ./ grid_spacing) ./ 2.355;
%     smoo = imgaussfilt3(M.T,sig,'padding','replicate');
%     MA.T_sm(:,:,zinds) = smoo(:,:,zinds);
% end
% 
% % ma = load('/Volumes/SDBlue/data/sgwex/3DST/Model_as_AIRS/20150705-1700_Model_as_AIRS_3DST_4dp.mat');















return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTY WOTTY PLOT PLOT - 2x2 panel with just raw Temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.625 1])

%-------------------------------------------------------
vert_gap    = 0.13;    horz_gap    = 0.105;
lower_marg  = 0.12;    upper_marg  = 0.025;
left_marg   = 0.085;   right_marg  = 0.055;

rows = 2; cols = 2;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 18;

%--------------------------------------------------------

drawmosaic = 1;
markersize = 8;

% define boundaries
% xbox = clon + [-6  8];
% ybox = clat + [-3.85 4];
xbox = [-42.1 -29.9];
ybox = [-58.1 -50.5];


TBP = {...
    A.Tg_lev-4,...
    M.T_lev,...
    MA.Tg_lev,...
    MA.Tpg_lev};

tbpxy = {...
    Airs.Lon,Airs.Lat;...
    Model.Lon,Model.Lat;...
    Airs.Lon,Airs.Lat;...
    Airs.Lon,Airs.Lat};

mosaics = {
    1,...
    0,...
    1,...
    1};

clims = {...
    243+plusminus(20),...
    245+plusminus(30),...
    242+plusminus(20),...
    [-15 15]};

cbarlabels = {...
    'Temperature (K)',...
    'Temperature (K)',...
    'Temperature (K)',...
    'Temperature (K)'};

cbarticks = {...
    200:10:300,...
    200:10:300,...
    200:10:300,...
    -20:5:20};

cbarlims = {...
    clims{1},...
    clims{2},...
    clims{3}+1,...
    clims{4}};

tits = {...
    'T',...
    'T''',...
    'T',...
    'T'''};



for ax = 1:3
    
    subplot(rows,cols,ax)
    
    axx = gca;
    
    % set limits first to speed up mosaic
    xlim(xbox)
    ylim(ybox)
    
    if drawmosaic
        switch mosaics{ax}
            case 1
                % MOSAIC PLOT
                data = TBP{ax};
                cl = clims{ax}; clim(cl); nclev = 31;
                cmap = cbrew('RdBu',nclev);
                Fcmap = griddedInterpolant({linspace(cl(1),cl(2),nclev),1:3},cmap,'linear','nearest');
                for i = 1:90
                    for j = 1:270
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
        % STANDARD PCOLOR
        hold on; pcolor(tbpxy{ax,1},tbpxy{ax,2},TBP{ax}); shat;
    end
    
    
    %%%% FORMATTING
    xtick(-50:2:-20);
    xminortick(-50:1:-20)
    ytick(-60:-40);
    yminortick(-60:0.5:-40);
    
%     switch ax
%         case {1,2,3}
%             ylabel('Latitude')
%     end
%     xlabel('Longitude')
    axx.XTickLabel = lonlabels(axx.XTick);
    axx.YTickLabel = latlabels(axx.YTick);
    
    
    switch ax
            
        case {1,2,3}
            
            % COLORBAR
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
            cbar.Ruler.MinorTickValues = 200:5:300;
            
            % reposition
            cbar.Position = [1 1 1 1] .* cbar.Position;
            drawnow
            apos = get(gca,'position');
            cbar.Position(3) = 0.6*apos(3);
            cbar.Position(1) = apos(1) + 0.2*apos(3);
            cbar.Position(2) = cbar.Position(2) - 0.085;
            cbar.Position(4) = 1.2*cbar.Position(4);
    end
    
    % COASTLINE:
    C = nph_draw_coastline([xbox(1) ybox(1) ; xbox(2) ybox(2)],0,'noplot');
    for d = 1:length(C)
        if length(C(d).Lon) > 1
            hold on; plot(C(d).Lon,C(d).Lat,'w','linewi',1.5)
            hold on; plot(C(d).Lon,C(d).Lat,'k','linewi',1)
        end
    end
    
    % draw some scan track lines
    switch ax
        case {1,2,3}
            lontrack = Airs.Lon(45,:);
            lattrack = Airs.Lat(45,:);
            linespec = {'linewi',2,'linest','--','color',rgbtrip(0.5)};
            hold on; plot(lontrack,lattrack,linespec{:});
    end
    
    % and some cross section lines
    switch ax
        case {1,2,3}
            xtlon = Airs.Lon(:,119);
            xtlat = Airs.Lat(:,119);
            missrng = inrange(1:90,[44 52]);
            linespec = {'linewi',1.5,'linest','--','color',[1 0 1 0.8]};
            % a full one missing out the main GWs by the island
            xtlon(missrng) = NaN;
            hold on; p = plot(xtlon,xtlat,linespec{:});
            % now a low alpha one for those perts
            xtlon = Airs.Lon(:,119);
            hold on; p = plot(xtlon(missrng),xtlat(missrng),linespec{:},'color',[1 0 1 0.3]);
    end
    
    
    
    clim(clims{ax})
    colormap(gca,cbrew('RdBu',31))
    
    % titles and letters
    drawnow;
    nph_text([0.025 0.87],['(' alphabet(ax) ')'],'fontsize',1.8*fs,'HorizontalAlignment','center','textborder','w');
    %     nph_text([0.1 0.865],tits{ax},'fontsize',1.6*fs,'fontweight','bold','HorizontalAlignment','center','textborder','w');
    
    setfont(fs)
    
    set(gca,'linewi',1.5,'color',rgbtrip(.95),'layer','top','tickdir','out')
    
    % outer line
    hold on; plot(xbox([1 2 2 1 1]),ybox([1 1 2 2 1]),'color','k','linewi',1.5);
    
end % next axes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% HORIZONTAL CUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xcut = 119;
A.Th  = A.Tg_lev(:,xcut);
MA.Th = MA.Tg_lev(:,xcut);

% the model is a little more difficult:
striplon = Airs.Lon(:,xcut); striplat = Airs.Lat(:,xcut);
% method: use a 5th order poly to fit the lat lon line of this cross track
% row.
p = polyfit(striplon,striplat,5);
% now evaluate this polynomial at super high longitude spacing. it won't be
% perfectly regularly spaced but should be very good:
striploni = linspace(Airs.Lon(1,xcut),Airs.Lon(90,xcut),2000);
striplati = polyval(p,striploni);
% now evaluate this through the full res model temperature:
% mt = movmean(M.T,9,1);
% F_m = scatteredInterpolant(Model.Lon(:),Model.Lat(:),linearise(mt(:,:,zind)),'linear','none');
M.Thi = F_m(striploni,striplati);
% also do some wind too:
M.windi = F_wind(striploni,striplati);
% find distance from nadir:
[d,~] = distance(Airs.Lat(45,xcut),Airs.Lon(45,xcut),striplat,striplon);
[di,~] = distance(Airs.Lat(45,xcut),Airs.Lon(45,xcut),striplati,striploni);
d = deg2km(d);
di = deg2km(di);
% set first half to be zero
dd = diff(d([1:end end])); ddi = diff(di([1:end end]));
d(dd < 0) = -d(dd < 0); di(ddi < 0) = -di(ddi < 0);

%%%% TOPOGRAPHY CUT
% trim to southern tip of the island: (done manually)
latrng = inrange(topolat,[-54.8 -54.64]);
topolat = topolat(latrng);
topo = double(topo(latrng,:));
topoline = nanmax(topo,[],1) ./ 1000; % convert to km
[dtp,~] = distance(Airs.Lat(45,xcut),Airs.Lon(45,xcut),Airs.Lat(45,xcut),topolon);
dtp = deg2km(dtp);
ddtp = diff(dtp([1:end end]));
dtp(ddtp < 0) = -dtp(ddtp < 0);




subplot(rows,cols,4)

axx = gca;
grid on

% Formatting first...
xlim([-200 250])
ylim([170 310])

xtick(-600:50:600);
xminortick(-600:25:600);

ytick(150:20:400);
yminortick(-150:10:400);

yl = ylabel('Temperature (K)');
yl.Position(1) = 1.1*yl.Position(1);

xlabel('Across-track distance (km)')

set(gca,'tickdir','out','linewi',1.5,'layer','top')

% nadir line
linespec = {'linewi',2,'linest','--','color',rgbtrip(0.5)};
hold on; plot([0 0],axx.YLim,linespec{:});

    % cut through line
    linespec = {'linewi',1.5,'linest','--','color',[1 0 1 0.8]};
    hold on; plot(axx.XLim,[240 240],linespec{:});

% outer line
hold on; plot(axx.XLim([1 2 2 1 1]),axx.YLim([1 1 2 2 1]),'color','k','linewi',1.5);

%%%% DATA
% CROSS SECTIONS
linespec = {'linewi',3};
airscolor = mean([mcolor(2) ; mcolor(7)]);
hold on; plot(di,M.Thi,linespec{:},'color','w','linewi',4);
hold on; M.plot  = plot(di,M.Thi,linespec{:},'color',[rgbtrip(.2) 0.6]);

hold on; plot(d,MA.Th,linespec{:},'color','w','linewi',4);
hold on; MA.plot = plot(d,MA.Th,linespec{:},'color',mcolor(1));

hold on; plot(d,A.Th-4,linespec{:},'color','w','linewi',4);
hold on; A.plot  = plot(d,A.Th-4,linespec{:},'color',airscolor);


% TOPOGRAPHY ON A NEW Y AXIS
yyaxis right;
hold on; fill([dtp reverse(dtp)],[topoline zeros(size(topoline))],'k');
ylim([0 70])
ytick([0 10]);
axx.YMinorTick = 'on';
axx.YAxis(2).MinorTickValues = 0:2.5:10;
set(gca,'ycolor','k','linewi',1.5)
ylabel('Altitude (km)                    ','HorizontalAlignment','right')

setfont(fs)


nph_text([0.01 0.87],['(' alphabet(4) ')'],'fontsize',1.8*fs,'HorizontalAlignment','center','textborder','w');


%%%% LEGEND
legend([M.plot MA.plot A.plot],{'Model','Model-as-AIRS','AIRS'})



return



%% EXPORT? ================================================================

savename = ['~/Desktop/' granule '_Mosaic_v2'];

disp(['Exporting to "' savename '"...'])
nph_saveas(gcf,savename,'png')
disp('Done.')






return






































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTY WOTTY PLOT PLOT - old version with T'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.675 1])

%-------------------------------------------------------
vert_gap    = 0.08;    horz_gap    = 0.12;
lower_marg  = 0.16;    upper_marg  = 0.025;
left_marg   = 0.085;   right_marg  = 0.05;

rows = 2; cols = 2;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 20;

%--------------------------------------------------------

drawmosaic = 0;
markersize = 8.3;

% define boundaries
% xbox = clon + [-6  8];
% ybox = clat + [-3.85 4];
xbox = [-42.1 -29.9];
ybox = [-58.1 -50.5];


TBP = {...
    A.Tg_lev,...
    A.Tpg_lev,...
    M.T_lev,...
    MA.Tpg_lev};

tbpxy = {...
    Airs.Lon,Airs.Lat;...
    Airs.Lon,Airs.Lat;...
    Model.Lon,Model.Lat;...
    Airs.Lon,Airs.Lat};

mosaics = {
    1,...
    1,...
    0,...
    1};

clims = {...
    [222.5 272.5],...
    [-15 15],...
    [210 280],...
    [-15 15]};

cbarlabels = {...
    'Temperature (K)',...
    'GW Perturbations (K)',...
    'Temperature (K)',...
    'GW Perturbations (K)'};

cbarticks = {...
    200:10:300,...
    -20:5:20,...
    200:10:300,...
    -20:5:20};

tits = {...
    'T',...
    'T''',...
    'T',...
    'T'''};



for ax = 1:1
    
    subplot(rows,cols,ax)
    
    % set limits first to speed up mosaic
    xlim(xbox)
    ylim(ybox)
    
    if drawmosaic
        switch mosaics{ax}
            case 1
                % MOSAIC PLOT
                data = TBP{ax};
                cl = clims{ax}; clim(cl); nclev = 31;
                cmap = cbrew('RdBu',nclev);
                Fcmap = griddedInterpolant({linspace(cl(1),cl(2),nclev),1:3},cmap,'linear','nearest');
                for i = 1:90
                    for j = 1:270
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
        % STANDARD PCOLOR
        hold on; pcolor(tbpxy{ax,1},tbpxy{ax,2},TBP{ax}); shat;
    end
    
    
    %%%% FORMATTING
    xtick(-50:2:-20);
    xminortick(-50:1:-20)
    ytick(-60:-40);
    yminortick(-60:0.5:-40);
    
    
    switch ax
        case {1,3}
            ylabel('Latitude')
    end
    
    
    
    switch ax
        case {1,2}
            
            xlabel('Longitude')
            
        case {3,4}
            
            % COLORBAR
            cbar = colorbar('southoutside');
            
            cbar.Label.String = cbarlabels{ax};
            cbar.Label.Rotation = 0;
            %             cbar.Label.VerticalAlignment = 'middle';
            %             cbar.Label.HorizontalAlignment = 'left';
            cbar.TickDirection = 'out';
            cbar.Ticks = cbarticks{ax};
            cbar.LineWidth = 1.5;
            
            % reposition
            cbar.Position = [1 1 1 1] .* cbar.Position;
            drawnow
            apos = get(gca,'position');
            cbar.Position(3) = 0.6*apos(3);
            cbar.Position(1) = apos(1) + 0.2*apos(3);
            cbar.Position(2) = cbar.Position(2) - 0.09;
            cbar.Position(4) = 1.2*cbar.Position(4);
    end
    
    % COASTLINE:
    C = nph_draw_coastline([xbox(1) ybox(1) ; xbox(2) ybox(2)],0,'noplot');
    for d = 1:length(C)
        if length(C(d).Lon) > 1
            hold on; plot(C(d).Lon,C(d).Lat,'w','linewi',1.5)
            hold on; plot(C(d).Lon,C(d).Lat,'k','linewi',1)
        end
    end
    
    % draw some scan track lines
    switch ax
        case {1,2}
            lontrack = Airs.Lon(45,:);
            lattrack = Airs.Lat(45,:);
            linespec = {'linewi',2,'linest','--','color',rgbtrip(0.5)};
            hold on; plot(lontrack,lattrack,linespec{:});
    end
    
    % and some cross section lines
    switch ax
        case {1,3}
            xtlon = Airs.Lon(:,119);
            xtlat = Airs.Lat(:,119);
            if ax == 2 || ax == 4
                xtlon(42:52) = NaN;
            end
            linespec = {'linewi',1.5,'linest','--','color',[1 0 1 0.8]};
            hold on; p = plot(xtlon,xtlat,linespec{:});
    end
    
    
    
    clim(clims{ax})
    colormap(gca,cbrew('RdBu',31))
    
    % titles and letters
    drawnow;
    nph_text([0.025 0.87],['(' alphabet(ax) ')'],'fontsize',1.8*fs,'HorizontalAlignment','center','textborder','w');
    %     nph_text([0.1 0.865],tits{ax},'fontsize',1.6*fs,'fontweight','bold','HorizontalAlignment','center','textborder','w');
    
    setfont(fs)
    
    set(gca,'linewi',1.5,'color',rgbtrip(.95),'layer','top','tickdir','out')
    
    % outer line
    hold on; plot(xbox([1 2 2 1 1]),ybox([1 1 2 2 1]),'color','k','linewi',1.5);
    
    
end % next axes



return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOW A LINE PLOT OF SOME CROSS-SECTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.675 1])

%-------------------------------------------------------
% vert_gap    = 0.15;    horz_gap    = 0.15;
% lower_marg  = 0.16;    upper_marg  = 0.025;
% left_marg   = 0.085;   right_marg  = 0.05;
% 
% rows = 2; cols = 3;
vert_gap    = 0.08;    horz_gap    = 0.12;
lower_marg  = 0.16;    upper_marg  = 0.025;
left_marg   = 0.085;   right_marg  = 0.05;

rows = 2; cols = 2;


subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 20;

%--------------------------------------------------------

axs = {[4 5],6};

%%%% GENERATE HORIZONTAL CUTS
xcut = 119;
A.Th  = A.Tg_lev(:,xcut);
MA.Th = MA.Tg_lev(:,xcut);

% the model is a little more difficult:
striplon = Airs.Lon(:,xcut); striplat = Airs.Lat(:,xcut);
% method: use a 5th order poly to fit the lat lon line of this cross track
% row.
p = polyfit(striplon,striplat,5);
% now evaluate this polynomial at super high longitude spacing. it won't be
% perfectly regularly spaced but should be very good:
striploni = linspace(Airs.Lon(1,xcut),Airs.Lon(90,xcut),2000);
striplati = polyval(p,striploni);
% now evaluate this through the full res model temperature:
% mt = movmean(M.T,9,1);
% F_m = scatteredInterpolant(Model.Lon(:),Model.Lat(:),linearise(mt(:,:,zind)),'linear','none');
M.Thi = F_m(striploni,striplati);
% % % % % also evaluate the topography too:
% % % % topo_xcut = F_topo(striploni,striplati);
% find distance from nadir:
[d,~] = distance(Airs.Lat(45,xcut),Airs.Lon(45,xcut),striplat,striplon);
[di,~] = distance(Airs.Lat(45,xcut),Airs.Lon(45,xcut),striplati,striploni);
d = deg2km(d);
di = deg2km(di);
% set first half to be zero
dd = diff(d([1:end end])); ddi = diff(di([1:end end]));
d(dd < 0) = -d(dd < 0); di(ddi < 0) = -di(ddi < 0);

%%%% TOPOGRAPHY CUT
% trim to southern tip of the island: (done manually)
latrng = inrange(topolat,[-54.8 -54.64]);
topolat = topolat(latrng);
topo = double(topo(latrng,:));
topoline = nanmax(topo,[],1) ./ 1000; % convert to km
[dtp,~] = distance(Airs.Lat(45,xcut),Airs.Lon(45,xcut),Airs.Lat(45,xcut),topolon);
dtp = deg2km(dtp);
ddtp = diff(dtp([1:end end]));
dtp(ddtp < 0) = -dtp(ddtp < 0);


%%%% GENERATE VERTICAL CUTS
airs_vert_loc = [48 xcut];
A.Tv = sq(A.Tg(airs_vert_loc(1),airs_vert_loc(2),:)); % that was easy

% for models, find distance from the airs nadir
[dv,~] = distance(Airs.Lat(airs_vert_loc(1),airs_vert_loc(2)),Airs.Lon(airs_vert_loc(1),airs_vert_loc(2)),Model.Lat,Model.Lon);
[~,ind] = min(abs(dv(:)));
[i,j] = ind2sub(size(Model.Lon),ind);
model_vert_loc = [i j];
nudge = [0 0];
% raw temps
M.Tv = sq(M.T(model_vert_loc(1)+nudge(1),model_vert_loc(2)+nudge(2),:));
MA.Tv = sq(MA.T_sm(model_vert_loc(1)+nudge(1),model_vert_loc(2)+nudge(2),:));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CROSS-SECTION PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% subplot(rows,cols,axs{1})
subplot(rows,cols,4)

axx = gca;
grid on

% Formatting first...
xlim([-250 250])
ylim([170 310])

xtick(-600:50:600);
xminortick(-600:25:600);

ytick(150:20:400);
yminortick(-150:10:400);

yl = ylabel('Temperature (K)');
yl.Position(1) = 1.1*yl.Position(1);

xlabel('Across-track distance (km)')

set(gca,'tickdir','out','linewi',1.5,'layer','top')

% nadir line
linespec = {'linewi',2,'linest','--','color',rgbtrip(0.5)};
hold on; plot([0 0],axx.YLim,linespec{:});

    % cut through line
    linespec = {'linewi',1.5,'linest','--','color',[1 0 1 0.8]};
    hold on; plot(axx.XLim,[240 240],linespec{:});

% outer line
hold on; plot(axx.XLim([1 2 2 1 1]),axx.YLim([1 1 2 2 1]),'color','k','linewi',1.5);

%%%% DATA
% CROSS SECTIONS
linespec = {'linewi',3};
airscolor = mean([mcolor(2) ; mcolor(7)]);
hold on; plot(di,M.Thi,linespec{:},'color','w','linewi',4);
hold on; M.plot  = plot(di,M.Thi,linespec{:},'color',[rgbtrip(.2) 0.6]);

hold on; plot(d,MA.Th,linespec{:},'color','w','linewi',4);
hold on; MA.plot = plot(d,MA.Th,linespec{:},'color',mcolor(1));

hold on; plot(d,A.Th-4,linespec{:},'color','w','linewi',4);
hold on; A.plot  = plot(d,A.Th-4,linespec{:},'color',airscolor);


% TOPOGRAPHY ON A NEW Y AXIS
yyaxis right;
hold on; fill([dtp reverse(dtp)],[topoline zeros(size(topoline))],'k');
ylim([0 70])
ytick([0 10]);
axx.YMinorTick = 'on';
axx.YAxis(2).MinorTickValues = 0:2.5:10;
set(gca,'ycolor','k','linewi',1.5)
ylabel('Altitude (km)             ','HorizontalAlignment','right')

setfont(fs)


nph_text([0.01 0.87],['(' alphabet(5) ')'],'fontsize',1.8*fs,'HorizontalAlignment','center','textborder','w');


%%%% LEGEND
legend([M.plot MA.plot A.plot],{'Model','Model-as-AIRS','AIRS'})



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% VERTICAL PROFILE PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

subplot(rows,cols,axs{2})

axx = gca;
grid on

% Formatting first...
xlim([170 310])
ylim([0 70])

xtick(-100:40:400);
xminortick(-100:10:400);

ytick(0:10:70);
yminortick(0:5:70);


xlabel('Temperature (K)')
ylabel({'Altitude (km)',' '})

set(gca,'tickdir','out','linewi',1.5,'layer','top')

% nadir line
linespec = {'linewi',2,'linest','--','color',rgbtrip(0.5)};
hold on; plot([0 0],axx.YLim,linespec{:});

% cut through line
linespec = {'linewi',1.5,'linest','--','color',[1 0 1 0.8]};
hold on; plot(axx.XLim,[45 45],linespec{:});

% outer line
hold on; plot(axx.XLim([1 2 2 1 1]),axx.YLim([1 1 2 2 1]),'color','k','linewi',1.5);

%%%% DATA
% CROSS SECTIONS
linespec = {'linewi',3};
airscolor = mean([mcolor(2) ; mcolor(7)]);

hold on; plot(M.Tv,Model.Alt,linespec{:},'color',[rgbtrip(.2) 0.6]);
hold on; plot(MA.Tv,Model.Alt,linespec{:},'color',mcolor(1));

atv = A.Tv-4; matv = nanmean(atv);
atv = (0.9 .* (atv - matv)) + matv;
hold on; plot(atv,Airs.Alt,linespec{:},'color',airscolor);

setfont(fs)





return



%% EXPORT? ================================================================

savename = ['~/Desktop/' granule '_Mosaic_pt2'];

disp(['Exporting to "' savename '"...'])
nph_saveas(gcf,savename,'png')
disp('Done.')




return
































































%%



files = dir([direc '*.mat']);

figure; hold all; whitefig; figpos([0.5 0.5])


% for i = 130:length(files)
for i = 131 % 20130710-0300
    
    
    granule = files(i).name(1:13);
    yrmonth = granule(1:6);
    
    if strcmpi(yrmonth,'201501')
        continue
    end
    
    cla
    
    disp(granule)
    
    load([direc files(i).name])
    
    [tp,bg] = nph_airs_4dp_detrend(Airs.T,1,4);
    
    zlev = 16;
    
    t_a = Airs.T(:,:,zlev);
    tp_a = tp(:,:,zlev);
    %     tp_lev = smoothn(tp_lev,[3 3]);
    
    hold on; pcolor(Airs.Lon,Airs.Lat,tp_a); shat;
    %     hold on; imagesc(Airs.Lon,Airs.Lat,flipud(tp_lev)); shat;
    %     hold on; pcolor(Airs.Lon,Airs.Lat,bg(:,:,zlev)); shat;
    
    % SGWEX MODEL CENTRE POINT:
    clat = -54.5269874182436;
    clon = -37.1367495437688;
    
    xbox = clon + [-7  8];
    ybox = clat + [-3.85 4];
    
    % COASTLINE:
    C = nph_draw_coastline([xbox(1) ybox(1) ; xbox(2) ybox(2)],0,'noplot');
    for d = 1:length(C)
        if length(C(d).Lon) > 10
            hold on; plot(C(d).Lon,C(d).Lat,'k','linewi',1.25)
        end
    end
    
    xlim(xbox)
    ylim(ybox)
    
    
    clim([-15 15])
    
    cbar = colorbar;
    colormap(gca,cbrew('RdBu',31))
    setfont(16)
    
    title(granule)
    
    drawnow
    
    %   pause
    
end



return



%% TRY AIRS AS A DOT MATRIX MAP

figure; hold all; whitefig; figpos([0.5 0.6])

%     hold on; pcolor(Airs.Lon,Airs.Lat,tp_lev); shat;

data = tp_a;

markersize = 11;

clims = [-15 15];
%     clims = [225 265];
clim(clims);
nclev = 31;

cmap = cbrew('RdBu',nclev);
Fcmap = griddedInterpolant({linspace(clims(1),clims(2),nclev),1:3},cmap,'linear','nearest');

xlim(xbox)
ylim(ybox)


for i = 1:90
    for j = 1:270
        % if in the region:
        if inrange(Airs.Lon(i,j),xbox) && inrange(Airs.Lat(i,j),ybox)
            if ~isnan(data(i,j))
                % get dot color:
                dotcolor = Fcmap({data(i,j),1:3});
                % draw a dot:
                hold on; plot(Airs.Lon(i,j),Airs.Lat(i,j),'marker','o','markerfacecolor',dotcolor,'markersize',markersize,'markeredgecolor','none');
            end
        end
    end
    drawnow;
end

% COASTLINE:
C = nph_draw_coastline([xbox(1) ybox(1) ; xbox(2) ybox(2)],0,'noplot');
for d = 1:length(C)
    if length(C(d).Lon) > 10
        hold on; plot(C(d).Lon,C(d).Lat,'k','linewi',1.25)
    end
end

cbar = nph_colorbar;
colormap(gca,cbrew('RdBu',31))
setfont(16)

set(gca,'color',rgbtrip(.9))



title(granule)



%%

% figure; hold all; whitefig;
%
% % tp_lev(43:52,119) = NaN
% rng = 44:52;
% loc = 119;
%
% % for i = 1:length(rng)
% %     hold on; plot(sq(tp(rng(i),loc,:))',Airs.Alt); shat
% %     pause
% % end
%
% % for z = 12:18
% %     if isodd(z)
% %         color = mcolor(1);
% %     else
% %         color = mcolor(2);
% %     end
% %     hold on; plot(1:90,Airs.T(:,loc,z),'color',color); grid on
% %     drawnow
% % %     pause
% % end
%
% hold on; plot(1:90,Airs.T(:,loc,17),'color',mcolor(1)); grid on
% hold on; plot(1:90,bg(:,loc,17),'color',mcolor(2)); grid on


% hold on; imagesc(rng,Airs.Alt,sq(tp(rng,loc,:))'); ydir
% hold on; contourf(rng,Airs.Alt,sq(tp(rng,loc,:))',-25:1:25,'edgecolor','none');


% clim([-20 20])
% cbar = colorbar;
%    colormap(gca,cbrew('RdBu',31))
%    setfont(16)



%
% hold on; pcolor(Airs.Lon,Airs.Lat,tp_lev); shat;
%
% xlim(xbox)
% ylim(ybox)
%
%
% clim([-15 15])
%
% cbar = colorbar;
% colormap(gca,cbrew('RdBu',31))
% setfont(16)





%% PUT MODEL ON AIRS GRID FOR THIS EXAMPLE?

load(['/Volumes/SDBlue/data/sgwex/3DST/Model/' granule '_Model_3DST_4dp.mat'])

%%

T = double(Model.T);

% simulate the airs footprint AND vertical resolution:
% MODEL AS AIRS ON AIRS GRID
% vert_res = airs3d_vert_res(Airs.Alt(zlev));
grid_spacing = [1.5 1.5 1.5];
sig = ([13.5 13.5 vert_res] ./ grid_spacing) ./ 2.355; % see hoffmann 2014, footprint size = 13.5km
T_ma = imgaussfilt3(T,sig,'padding','replicate');

% MODEL ON AIRS GRID
grid_spacing = [1.5 1.5];
sig = ([13.5 13.5] ./ grid_spacing) ./ 2.355; % see hoffmann 2014, footprint size = 13.5km
T_m = imgaussfilt(T,sig,'padding','replicate');

% extract level
[~,zind] = min((Model.Alt - Airs.Alt(zlev)).^2);

% put on airs grid:
F_ma = scatteredInterpolant(Model.Lon(:),Model.Lat(:),linearise(T_ma(:,:,zind)),'linear','none');
t_ma = F_ma(Airs.Lon,Airs.Lat);

F_m = scatteredInterpolant(Model.Lon(:),Model.Lat(:),linearise(T_m(:,:,zind)),'linear','none');
t_m = F_m(Airs.Lon,Airs.Lat);

% get perts on cross track grid
[tp_ma,bg_ma] = nph_airs_4dp_detrend(t_ma,1,4);
[tp_m,bg_m] = nph_airs_4dp_detrend(t_m,1,4);


%% MODEL ON AIRS GRID

figure; hold all; whitefig; figpos([0.5 0.5])

hold on; pcolor(Airs.Lon,Airs.Lat,tp_m); shat;


% COASTLINE:
C = nph_draw_coastline([xbox(1) ybox(1) ; xbox(2) ybox(2)],0,'noplot');
for d = 1:length(C)
    if length(C(d).Lon) > 10
        hold on; plot(C(d).Lon,C(d).Lat,'k','linewi',1.25)
    end
end

xlim(xbox)
ylim(ybox)

clim([-15 15])

cbar = colorbar;
colormap(gca,cbrew('RdBu',31))
setfont(16)

title('Model on AIRS grid')


%% MODEL AS AIRS ON AIRS GRID

figure; hold all; whitefig; figpos([0.5 0.5])

hold on; pcolor(Airs.Lon,Airs.Lat,tp_ma); shat;

% COASTLINE:
C = nph_draw_coastline([xbox(1) ybox(1) ; xbox(2) ybox(2)],0,'noplot');
for d = 1:length(C)
    if length(C(d).Lon) > 10
        hold on; plot(C(d).Lon,C(d).Lat,'k','linewi',1.25)
    end
end

xlim(xbox)
ylim(ybox)

clim([-15 15])

cbar = colorbar;
colormap(gca,cbrew('RdBu',31))
setfont(16)

title('Model as AIRS on AIRS grid')


%% across-track temps:

figure; hold all; whitefig;

% tp_lev(43:52,119) = NaN
% rng = 44:52;
% loc = 119;

% hold on; plot(1:90,Airs.T(:,119,17),'color',mcolor(1)); grid on
% hold on; plot(1:90,t_mod(:,118),'color',mcolor(2)); grid on

linespec = {'linewi',2};

hold on; plot(1:90,M.Tpg_lev(:,119),'color',[rgbtrip(.2) 0.5],linespec{:}); grid on
hold on; plot(1:90,A.Tpg_lev(:,119),'color',mcolor(2),linespec{:}); grid on
hold on; plot(1:90,MA.Tpg_lev(:,119),'color',mcolor(1),linespec{:}); grid on

xlim([10 80])


%% full res model pics

figure; pcolor(Model.Lon,Model.Lat,Model.T(:,:,zind)); shat

clim([200 290])

cbar = colorbar;
colormap(gca,cbrew('RdBu',31))
setfont(16)

% COASTLINE:
C = nph_draw_coastline([xbox(1) ybox(1) ; xbox(2) ybox(2)],0,'noplot');
for d = 1:length(C)
    if length(C(d).Lon) > 1
        hold on; plot(C(d).Lon,C(d).Lat,'k','linewi',1.25)
    end
end

xlim(xbox)
ylim(ybox)



































