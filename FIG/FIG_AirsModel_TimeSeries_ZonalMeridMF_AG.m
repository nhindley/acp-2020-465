



% Aaaaaand and even newer version of this for the airs-gridded (AG) version
% of the model-as-AIRS.




% - Even Newer Timeseries of AIRS/Model Amplitudes and Zonal and Meridional
%   MF against time, with directional wind speeds below.
%   This version does not do different heights, but does do median and
%   percentiles for wave amplitudes.
%   This version probably superceeds the other timeseries programs now.

% - Also the distribution of wave amplitudes with horizontal distance from
%   the island.




runagain = 1;


%% DIRECTORIES ============================================================

% matlabdirec         = '/Users/neil/Drive/MATLAB/';
% airsdirec           = '/Volumes/SDBlue/data/sgwex/3DST/Airs/downsized_new/';
% modeldirec          = '/Volumes/SDBlue/data/sgwex/3DST/Model/downsized_new/';
% modelasairsdirec    = '/Volumes/SDBlue/data/sgwex/3DST/Model_as_AIRS/downsized_new/';
% winddirec           = '/Volumes/SDBlue/data/sgwex/VertWindProfiles/';
%
% matlabdirec         = '/Users/neil/Drive/MATLAB/';
% airsdirec           = '/Volumes/SDBlue/data/sgwex/3DST/Airs/4dp/downsized/';
% modeldirec          = '/Volumes/SDBlue/data/sgwex/3DST/Model/4dp/downsized/';
% modelasairsdirec    = '/Volumes/SDBlue/data/sgwex/3DST/Model_as_AIRS/4dp/downsized/';
% winddirec           = '/Volumes/SDBlue/data/sgwex/VertWindProfiles/';

matlabdirec         = '/Users/neil/Drive/MATLAB/';
airsdirec           = '/Volumes/SDBlue/data/sgwex/3DST/Airs/AG_sm/';
modelasairsdirec    = '/Volumes/SDBlue/data/sgwex/3DST/Model_as_AIRS/AG_sm/';
winddirec           = '/Volumes/SDBlue/data/sgwex/VertWindProfiles/';


% %% LOAD TOPO AND MAP ======================================================
% % need this for the amp/gini/distance plots :)
%
% disp('Loading Topography...')
% load('/Users/neil/Drive/MATLAB/topography/SRTM/srtm_29_23/srtm_29_23.mat')
%
% topo = single(Topo.Topo);
% topo(topo < 0) = 0;
% topo = double(nanmax(topo,[],1));
%
clat = -54.5269874182436;
clon = -37.1367495437688;
%
% [d,az] = distance(clat,clon,repmat(clat,1,length(Topo.Lon)),Topo.Lon);
% xdist = deg2km(d).*sind(az);
%
% % smooooth
% topo = movmean(topo,51);

%% CREATE GRIDS AND REGIONS ===============================================

gsz = [90 120 50];
% gsz = [60 80 50];

vec1 = linspace(-450,450,gsz(1));
vec2 = linspace(-600,600,gsz(2));
vec3 = 1.5:1.5:75;
[D1,D2,D3] = ndgrid(vec1,vec2,vec3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DEFINE REGIONS

%                d1     r     d2
%  --------------.------------.--------
% |              |       '    |        |
% |                   '                |
% |      A       |  '       B | r      |
% |                '                   |
% |              |'           |        |
% | - - - - - - - - - SG- - - - - - -  | h
% |              |'           |        |
% |                '                   |
% |              |  '         | r      |
% |                   '                |
% |              |       '    |        |
%  ---------------------------'---------
%                    L

% Region A and B have equal area. Region B is the intersection of the
% circle of radius r centred at d1+r, and the rectangluar region h*d2.
% Region A is the exclusion of this.

% based on some maths I did on paper, d1 is this:

% model domain:
h = 400;
L = 600;

r = h/2;
d1 = (pi*h/8) - (h/4) + ((L-r)/2);
d2 = L - r - d1;

ring = quadadd(D1,D2-(d1+r-(L/2))) <= r;

zlims = [25 45];

% keep awayyy from the edges
edgebox = inrange(D1,pm(h/2)) & inrange(D2,pm(L/2));

% Region A - Upwind of the island
Regions.A = ~ring & edgebox & inrange(D2,[-(L/2) (L/2)-d2]) & inrange(D3,zlims);

% Region B - Over/Downwind of the island
Regions.B = (ring | inrange(D2,[(L/2)-d2 (L/2)])) & edgebox & inrange(D3,zlims);

% Region C - Both!
Regions.C = Regions.A | Regions.B;

R = {'A','B','C'};

% PERCENTILES
percentiles = [5 10 15 25 75 85 90 95];
for c = 1:length(percentiles)
    %     props{end+1} = ['p' num2str(percentiles(c))];
end

%% LOAD DATA ==============================================================

if runagain
    
    airsfiles = [...
        dir([airsdirec '201306*Airs*.mat']) ; ...
        dir([airsdirec '201307*Airs*.mat']) ; ...
        dir([airsdirec '201506*Airs*.mat']) ; ...
        dir([airsdirec '201507*Airs*.mat'])];
    
    modelasairsfiles = [...
        dir([modelasairsdirec '201306*Model*.mat']) ; ...
        dir([modelasairsdirec '201307*Model*.mat']) ; ...
        dir([modelasairsdirec '201506*Model*.mat']) ; ...
        dir([modelasairsdirec '201507*Model*.mat'])];
    
    files = {airsfiles,modelasairsfiles};
    direcs = {airsdirec,modelasairsdirec};
    
    
    % for each dataset...
    for q = 1:2
        
        files{q}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% ASSIGN OUTPUT STRUCTURE
        OUT = struct;
        
        for r = 1:3
            
            % Time
            OUT.Time = nan(1,length(files{q}));
            
            % Wave amplitude
            OUT.(R{r}).amedian = nan(1,length(files{q}));
            OUT.(R{r}).amean   = nan(1,length(files{q}));
            
            % Amplitude percentiles:
            for c = 1:length(percentiles)
                OUT.(R{r}).(['p' num2str(percentiles(c))]) = nan(1,length(files{q}));
            end
            
            % directional MF
            OUT.(R{r}).mf1    = nan(1,length(files{q}));
            OUT.(R{r}).mf1pos = nan(1,length(files{q}));
            OUT.(R{r}).mf1neg = nan(1,length(files{q}));
            
            OUT.(R{r}).mf2    = nan(1,length(files{q}));
            OUT.(R{r}).mf2pos = nan(1,length(files{q}));
            OUT.(R{r}).mf2neg = nan(1,length(files{q}));
            
            OUT.(R{r}).mf     = nan(1,length(files{q}));
            
            % MF percentiles
            for c = 1:length(percentiles)
                OUT.(R{r}).(['mf1p' num2str(percentiles(c))]) = nan(1,length(files{q}));
                OUT.(R{r}).(['mf2p' num2str(percentiles(c))]) = nan(1,length(files{q}));
            end
            
            % total MF
            OUT.(R{r}).totmf   = nan(1,length(files{q}));
            OUT.(R{r}).totmf_n = nan(1,length(files{q}));
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% for each timestep or overpass...
        for i = 1:length(files{q})
            
            disp(files{q}(i).name)
            
            granule = files{q}(i).name(1:13);
            
            %%%% load data...
            try
                load([direcs{q} files{q}(i).name])
            catch err
                err
                continue
            end
            
            % and assign:
            switch q
                case 1
                    Q = Airs;
                    %                     Q.Time = datenum(Q.Time);
                    Q.Time = Q.StartTime;
                    
                case 2
                    Q = Model;
                    Q.Time = datenum(Q.Time);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % apply UPWARD propagation assumption:
            posinds = Q.ST.F3 > 0;
            Q.ST.F1(posinds) = -Q.ST.F1(posinds);
            Q.ST.F2(posinds) = -Q.ST.F2(posinds);
            Q.ST.F3          = -abs(Q.ST.F3);
            Q.MF1(posinds)   = -Q.MF1(posinds);
            Q.MF2(posinds)   = -Q.MF2(posinds);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % now compute the data!!
            for r = 1:3
                
                inds = Regions.(R{r});
                
                % Extract Time
                OUT.Time(i) = Q.Time;
                
                % Wave amplitude
                amp = linearise(Q.ST.HA(inds));
                OUT.(R{r}).amedian(i) = nanmedian(amp);
                OUT.(R{r}).amean(i) = nanmean(amp);
                
                % amplitude percentiles:
                for c = 1:length(percentiles)
                    OUT.(R{r}).(['p' num2str(percentiles(c))])(i) = prctile(amp,percentiles(c));
                end
                
                % directional MF
                mf1pos = Q.MF1; mf1pos(Q.MF1 < 0) = 0;
                mf1neg = Q.MF1; mf1neg(Q.MF1 > 0) = 0;
                mf2pos = Q.MF2; mf2pos(Q.MF2 < 0) = 0;
                mf2neg = Q.MF2; mf2neg(Q.MF2 > 0) = 0;
                
                OUT.(R{r}).mf1pos(i) = nanmean(linearise(mf1pos(inds)));
                OUT.(R{r}).mf1neg(i) = nanmean(linearise(mf1neg(inds)));
                OUT.(R{r}).mf2pos(i) = nanmean(linearise(mf2pos(inds)));
                OUT.(R{r}).mf2neg(i) = nanmean(linearise(mf2neg(inds)));
                
                OUT.(R{r}).mf1(i)    = nanmean(linearise(Q.MF1(inds)));
                OUT.(R{r}).mf2(i)    = nanmean(linearise(Q.MF2(inds)));
                
                % abs MF
                mf = quadadd(Q.MF1,Q.MF2);
                OUT.(R{r}).mf(i)    = nanmean(linearise(mf(inds)));
                
                % total MF
                OUT.(R{r}).totmf(i)   = nansum(linearise(mf(inds)));
                OUT.(R{r}).totmf_n(i) = sum(double(~isnan(linearise(mf(inds)))));
                
                % MF percentiles:
                for c = 1:length(percentiles)
                    OUT.(R{r}).(['mf1p' num2str(percentiles(c))])(i) = prctile(Q.MF1(inds),percentiles(c));
                    OUT.(R{r}).(['mf2p' num2str(percentiles(c))])(i) = prctile(Q.MF2(inds),percentiles(c));
                end
                
            end
            
        end % next timestep or overpass...
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % finally, assign output:
        switch q
            case 1
                A = OUT;
            case 2
                MA = OUT;
        end
        
    end % next dataset
    
else
    
    %%%% LOAD
    load('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_TimeSeries_ZonalMeridMF_AG')
    disp('Loading saved data...')
    
    return
    
    %% SAVE?
    save('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_TimeSeries_ZonalMeridMF_AG','A','MA','Regions','h','L')
    disp('Saved.')
    
    return
    
    
end % end if runagain


return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PLOT TIMESERIES OF WAVE AMPLITUDE, ZONAL GWMF, MERIDIONAL GWMF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.9 1])

%-------------------------------------------------------
vert_gap    = 0.06;    horz_gap    = 0.06;
lower_marg  = 0.06;    upper_marg  = 0.05;
left_marg   = 0.085;     right_marg  = 0.065;

rows = 3; cols = 2;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 20;

%--------------------------------------------------------

xlims = {...
    [datenum(2013,07,01,16,00,00) datenum(2013,07,31,17,00,00)],...
    [datenum(2015,06,12,16,00,00) datenum(2015,07,07,03,00,00)]};

years = {2013,2015,2013,2015,2013,2015};

tits = {'AIRS','AIRS',' Model as AIRS',' Model as AIRS',' Model',' Model'};

%%%% PLOTTING SETTINGS
ylims = {...
    [0 6],[0 6],...
    [-80 20],[-80 20],...
    [-50 20],[-50 20],...
    };
ytix = {0:1:10,0:1:10,...
    -100:20:100,-100:20:100,...
    -100:10:100,-100:10:100,...
    };
yminortix = {0:0.5:10,0:0.5:10,...
    -100:10:100,-100:10:100,...
    -100:5:100,-100:5:100,...
    };
% ylabs = {'Wave Amplitude (K)',...
%     '{\bf{MF_x}} (mPa)',...
%     '{\bf{MF_y}} (mPa)',...
%     };
ylabs = {{'Wave Amplitude','(K)'},...
    {'Zonal MF','(mPa)'},...
    {'Meridional MF','(mPa)'},...
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tbp = {'','','mf2','mf2','mf1','mf1'};

acolor = [183 23 21]./255;
macolor = [33 114 180]./255;

for ax = 1:6
    
    yr = years{ax};
    yrstr = num2str(yr);
    
    subplot(rows,cols,ax)
    
    region = 'C';
    
    % smoothing?
    asmoo = 1;
    masmoo = 3;
    
    % LINESPEC
    msize = 24;
    airslinespec = {'linewi',2.5,'color',acolor,'marker','.','markersize',msize};
    modelasairslinespec = {'linewi',2.5,'color',macolor};
    
    % for each dataset...
    for q = 1:2
        switch q
            case 1
                Q = A;
                smoo = 1;
                edgecolor = 'w';
                facecolor = rgbtrip(.85);
                facealpha = 1;
            case 2
                Q = MA;
                smoo = 5;
                edgecolor = 'w';
                facecolor = rgbtrip(.6);
                facealpha = 1;
        end
        
        % filled regions
        ti = inrange(Q.Time,xlims{2-double(isodd(ax))} + [0 0]);
        tim = Q.Time(ti);
        
        switch ax
            case {1,2}
                pcs = [10 90];
            case {3,4,5,6}
                pcs = [10 90];
        end
        
        filldata1 = movmean(Q.(region).([tbp{ax} 'p' num2str(pcs(1))])(ti),smoo,2,'omitnan');
        filldata2 = movmean(Q.(region).([tbp{ax} 'p' num2str(pcs(2))])(ti),smoo,2,'omitnan');
        fillline = [filldata1 reverse(filldata2)];
        
        % FILLED AREAS
        switch ax
            case {1,2}
                hold on; F = fill([tim reverse(tim)],fillline,'k');
            case {3,4,5,6}
                hold on; F = fill([tim reverse(tim)],1000.*fillline,'k');
        end
        F.FaceColor = facecolor;
        F.FaceAlpha = facealpha;
        F.EdgeColor = edgecolor;
        F.EdgeAlpha = 0.75;
        F.LineWidth = 1;
        F.LineStyle = '-';
        
    end
    
    % LINE PLOTS
    for q = [1 2]
        
        switch q
            case 1
                Q = A;
                smoo = 1;
                linespec = airslinespec;
            case 2
                Q = MA;
                smoo = 3;
                linespec = modelasairslinespec;
        end
        
        % filled regions
        ti = inrange(Q.Time,xlims{2-double(isodd(ax))} + [0 0]);
        tim = Q.Time(ti);
        
        switch ax
            case {1,2}
                data = movmean(Q.(region).amean(ti),smoo,2,'omitnan');
            case {3,4,5,6}
                data = movmean(1000.*Q.(region).(tbp{ax})(ti),smoo,2,'omitnan');
        end
        
        switch q
            case 1 % A
                hold on; plot(tim,data,linespec{:},'color','w','linewi',3,'markersize',msize+4); grid on;
                hold on; plot(tim,data,linespec{:}); grid on;
            case 2 % MA
                hold on; plot(tim,data,linespec{:},'color','w','linewi',3.25); grid on;
                hold on; plot(tim,data,linespec{:}); grid on;
        end
    end
    
    
    % Axes formatting
    axx = gca;
    
    % Dates...
    xlim(xlims{2-mod(ax,2)})
    axx.XTick = [...
        datenum(yr,06,[1 5:5:30]) ...
        datenum(yr,07,[  5:5:30]) ...
        datenum(yr,08,[  5:5:25])];
    axx.XMinorTick = 'on';
    axx.XAxis.MinorTickValues = datenum(yr,06,1:100);
    datetick('x','dd mmmm','keepticks','keeplimits')
    
    ylim(ylims{ax})
    axx.YTick = ytix{ax};
    axx.YMinorTick = 'on';
    axx.YAxis.MinorTickValues = yminortix{ax};
    
    % YEAR LABEL
%     xlabel(yrstr,'fontweight','bold')
    switch ax
        case 1
            title('July 2013','fontweight','bold','fontsize',fs)
        case 2
            title('June - July 2015','fontweight','bold','fontsize',fs)    
    end
    
    set(gca,'tickdir','out','linewi',1.5,'clipping','off','layer','top')
    
    % YLABELS
    switch ax
        case {1,3,5}
            yl = ylabel(ylabs{ceil(ax/2)},'fontweight','bold');
            yl.Position(1) = yl.Position(1) - 1;
        case {2,4,6}
            yl = ylabel({ylabs{ceil(ax/2)}{2},' '},'fontweight','bold');
            yl.Position(1) = yl.Position(1) + 0.15;
    end
    
    drawnow;
%     nph_text([-0.02 1.05],['(' alphabet(ax) ')'],'fontsize',1.75*fs,'HorizontalAlignment','center','textborder','w');
    nph_text([0.01 0.8725],['(' alphabet(ax) ')'],'fontsize',1.7*fs,'HorizontalAlignment','left','textborder','w');
    drawnow;
    
    setfont(fs)
    
    %     return
    
end




return



%% EXPORT? ================================================================


savename = ['~/Desktop/AirsModel_Timeseries_A_MFx_MFy_AG_sm'];
disp(['Exporting to ' savename '...'])

nph_saveas(gcf,savename,'png')

disp('Done.')







return


%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TIMESERIES OF % OF GWMF UPWIND/DOWNWIND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.9 1])

%-------------------------------------------------------
vert_gap    = 0.06;    horz_gap    = 0.06;
lower_marg  = 0.06;    upper_marg  = 0.05;
left_marg   = 0.085;     right_marg  = 0.065;

rows = 3; cols = 2;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 20;

%--------------------------------------------------------

xlims = {...
    [datenum(2013,07,01,16,00,00) datenum(2013,07,31,17,00,00)],...
    [datenum(2015,06,12,16,00,00) datenum(2015,07,07,03,00,00)]};

years = {2013,2015,2013,2015,2013,2015};

tits = {'AIRS','AIRS',' Model as AIRS',' Model as AIRS',' Model',' Model'};

%%%% PLOTTING SETTINGS
% ylims = {...
%     [0.01 10],[0.01 10],...
%     [-80 20],[-80 20],...
%     [-50 20],[-50 20],...
%     };
% ytix = {10.^(-2:1:2),10.^(-2:1:2),...
%     -100:10:100,-100:10:100,...
%     -100:10:100,-100:10:100,...
%     };
% yminortix = {unique([0.01:0.01:0.1 0.1:0.1:1 1:10]),unique([0.01:0.01:0.1 0.1:0.1:1 1:10]),...
%     -100:5:100,-100:5:100,...
%     -100:5:100,-100:5:100,...
%     };
% ylabs = {{'Ratio of GWMF','A / B'},...
%     {'Zonal MF','[mPa]'},...
%     {'Meridional MF','[mPa]'},...
%     };

ylims = {...
    [0 100],[0 100],...
    [-80 20],[-80 20],...
    [-50 20],[-50 20],...
    };
ytix = {0:10:100,0:10:100,...
    -100:10:100,-100:10:100,...
    -100:10:100,-100:10:100,...
    };
yminortix = {0:10:100,0:10:100,...
    -100:5:100,-100:5:100,...
    -100:5:100,-100:5:100,...
    };
ylabs = {{'Fraction of MF','in Region B (%)'},...
    {'Zonal MF','[mPa]'},...
    {'Meridional MF','[mPa]'},...
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tbp = {'','','mf2','mf2','mf1','mf1'};

acolor = [183 23 21]./255;
macolor = [33 114 180]./255;

for ax = 1:2
    
    yr = years{ax};
    yrstr = num2str(yr);
    
    subplot(rows,cols,ax)
    
    region = 'C';
    
    % smoothing?
    asmoo = 1;
    masmoo = 3;
    
    % LINESPEC
    msize = 24;
    airslinespec = {'linewi',2.5,'color',acolor,'marker','.','markersize',msize};
    modelasairslinespec = {'linewi',2.5,'color',macolor};
    
    %%%% 'ZERO' LINE
    zeroline_x = min(xlims{ax}):max(xlims{ax});
    hold on; plot(zeroline_x,50*ones(size(zeroline_x)),'color',rgbtrip(0.4),'linewi',2,'linest','--'); 
    
    % LINE PLOTS
    for q = [2 1]
        
        switch q
            case 1
                Q = A;
                smoo = 1;
                linespec = airslinespec;
            case 2
                Q = MA;
                smoo = 5;
                linespec = modelasairslinespec;
        end
        
        % filled regions
        ti = inrange(Q.Time,xlims{2-double(isodd(ax))} + [0 0]);
        tim = Q.Time(ti);
        
        % compute the modifier to cope with missing data in each region:
        modifier = Q.A.totmf_n./Q.B.totmf_n;
        
        % frac of total in region B:
        frac = 100 .* (Q.B.totmf ./ Q.C.totmf) ./ modifier;
%         frac = (Q.A.totmf ./ Q.B.totmf) ./ modifier;
        
        switch q
            case 1 % A (also draw overpass times)
                hold on; plot(tim,frac(ti),linespec{:},'color','w','linewi',3,'markersize',msize+4); grid on;
                hold on; plot(tim,frac(ti),linespec{:}); grid on;
            case 2 % MA
                frac = movmean(frac,smoo,2,'omitnan');
                hold on; plot(tim,frac(ti),linespec{:},'color','w','linewi',3); grid on;
                hold on; plot(tim,frac(ti),linespec{:}); grid on;
        end
    end
    
%     logscale('y')
    
    % Axes formatting
    axx = gca;
    
    % Dates...
    xlim(xlims{2-mod(ax,2)})
    axx.XTick = [...
        datenum(yr,06,[1 5:5:30]) ...
        datenum(yr,07,[  5:5:30]) ...
        datenum(yr,08,[  5:5:25])];
    axx.XMinorTick = 'on';
    axx.XAxis.MinorTickValues = datenum(yr,06,1:100);
    datetick('x','dd mmmm','keepticks','keeplimits')
    
    ylim(ylims{ax})
    axx.YTick = ytix{ax};
    axx.YMinorTick = 'on';
    axx.YAxis.MinorTickValues = yminortix{ax};
    
    axx.YMinorGrid = 'on';
    
    % adjust yticklabels for log
    axx.YTickLabel = cellstr(num2str(axx.YTick'));
    
    % YEAR LABEL
%     xlabel(yrstr,'fontweight','bold')
%     switch ax
%         case 1
%             title('July 2013','fontweight','bold','fontsize',fs)
%         case 2
%             title('June - July 2015','fontweight','bold','fontsize',fs)    
%     end
    
    set(gca,'tickdir','out','linewi',1.5,'clipping','off','layer','top')
    
    % YLABELS
    switch ax
        case {1,3,5}
            yl = ylabel(ylabs{ceil(ax/2)},'fontweight','bold');
%             yl.Position(1) = yl.Position(1);
        case {2,4,6}
            yl = ylabel({'%',' '},'fontweight','bold');
            yl.Position(1) = yl.Position(1) + 0.15;
    end
    
    
    drawnow;
%     nph_text([-0.02 1.05],['(' alphabet(ax) ')'],'fontsize',1.75*fs,'HorizontalAlignment','center','textborder','w');
    nph_text([0.01 0.8725],['(' alphabet(ax+6) ')'],'fontsize',1.7*fs,'HorizontalAlignment','left','textborder','w');
    drawnow;
    
    setfont(fs)
    
    %     return
    
end




return



%% EXPORT? ================================================================


savename = ['~/Desktop/AirsModel_Timeseries_FracAB_v2'];
disp(['Exporting to ' savename '...'])

nph_saveas(gcf,savename,'png')

disp('Done.')







return


%%






































