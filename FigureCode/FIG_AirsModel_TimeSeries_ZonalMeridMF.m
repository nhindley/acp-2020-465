


% - Even Newer Timeseries of AIRS/Model Amplitudes and Zonal and Meridional
%   MF against time, with directional wind speeds below.
%   This version does not do different heights, but does do median and
%   percentiles for wave amplitudes.
%   This version probably superceeds the other timeseries programs now.

% - Also the distribution of wave amplitudes with horizontal distance from
%   the island.




runagain = 0;

load('/Users/neil/data/20150705-1700_AirsSG_3DST.mat')
sfvec = Airs.AltAmpScaling;
sf = permute(repmat(sfvec,1,60,80),[2 3 1]);



%% DIRECTORIES ============================================================

% matlabdirec         = '/Users/neil/Drive/MATLAB/';
% airsdirec           = '/Volumes/SDBlue/data/sgwex/3DST/Airs/downsized_new/';
% modeldirec          = '/Volumes/SDBlue/data/sgwex/3DST/Model/downsized_new/';
% modelasairsdirec    = '/Volumes/SDBlue/data/sgwex/3DST/Model_as_AIRS/downsized_new/';
% winddirec           = '/Volumes/SDBlue/data/sgwex/VertWindProfiles/';

matlabdirec         = '/Users/neil/Drive/MATLAB/';
airsdirec           = '/Volumes/SDBlue/data/sgwex/3DST/Airs/4dp/downsized/';
modeldirec          = '/Volumes/SDBlue/data/sgwex/3DST/Model/4dp/downsized/';
modelasairsdirec    = '/Volumes/SDBlue/data/sgwex/3DST/Model_as_AIRS/4dp/downsized/';
winddirec           = '/Volumes/SDBlue/data/sgwex/VertWindProfiles/';

%% LOAD TOPO AND MAP ======================================================
% need this for the amp/gini/distance plots :)

disp('Loading Topography...')
load('/Users/neil/Drive/MATLAB/topography/SRTM/srtm_29_23/srtm_29_23.mat')

topo = single(Topo.Topo);
topo(topo < 0) = 0;
topo = double(nanmax(topo,[],1));

clat = -54.5269874182436;
clon = -37.1367495437688;

[d,az] = distance(clat,clon,repmat(clat,1,length(Topo.Lon)),Topo.Lon);
xdist = deg2km(d).*sind(az);

% smooooth
topo = movmean(topo,51);



%% CREATE GRIDS ===========================================================

[D1,D2,D3] = ndgrid(...
    linspace(-450,450-15,60),...
    linspace(-600,600-15,80),...
    1.5:1.5:75);

d1 = linspace(-450,450-15,60);

% R = nanmean(quadadd(D1,D2),3);
% DD2 = nanmean(D2,3);

% define horz ring around the island:
r = 250; %km
xpos = 100; % km x-horizontally away from model centre

% height limits
zlims = [25 45];

% Get indexes for each height:
inds = quadadd(D1,D2-xpos) <= r & inrange(D3,zlims);
% inds2 = quadadd(D1,D2-xpos) <= 250 & inrange(D3,zlims);

% Distance from the island?
rbins = -800:25:800;
[XI,YI] = meshgrid(rbins,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if runagain
    
    airsfiles = [...
        dir([airsdirec '201306*AirsSG*.mat']) ; ...
        dir([airsdirec '201307*AirsSG*.mat']) ; ...
        dir([airsdirec '201506*AirsSG*.mat']) ; ...
        dir([airsdirec '201507*AirsSG*.mat'])];
    
    modelfiles = [...
        dir([modeldirec '201306*Model*.mat']) ; ...
        dir([modeldirec '201307*Model*.mat']) ; ...
        dir([modeldirec '201506*Model*.mat']) ; ...
        dir([modeldirec '201507*Model*.mat'])];
    
    modelasairsfiles = [...
        dir([modelasairsdirec '201306*Model*.mat']) ; ...
        dir([modelasairsdirec '201307*Model*.mat']) ; ...
        dir([modelasairsdirec '201506*Model*.mat']) ; ...
        dir([modelasairsdirec '201507*Model*.mat'])];
    
    files = {airsfiles,modelfiles,modelasairsfiles};
    direcs = {airsdirec,modeldirec,modelasairsdirec};
    
%     % bin all flux into horz wavelength bins
%     bins.lh = [0:25:250 300 400 600 1000];
%     [XI.lh,YI.lh] = meshgrid(bins.lh,1);
    
%     %and amplitude bins
%     bins.a = 0:0.25:25;
%     [XI.a,YI.a] = meshgrid(bins.a,1);
    
%     % props = {'lh','a','mf','mf1','mf1pos','mf1neg','mf1pos_median','mf1neg_median','fracmf1pos','fracmf1neg','mf2','Time'};
    props = {...
        'mf','mf1','mf2',...
        'mf1pos','mf1neg',...
        'mf2pos','mf2neg',...
        'a','astd','amedian',...
        'ad','mfd',...
        'Time'};
    
    percentiles = [5:5:95];
    for c = 1:length(percentiles)
        props{end+1} = ['p' num2str(percentiles(c))];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% LOAD AIRSO, MODEL AND MODELO AS AIRSO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % for each dataset...
    for q = 1:3
        
        % ASSIGN OUTPUT STRUCTURE
        OUT = struct;
        
        % less code to assign fields to structs in a loop:
        for p = 1:length(props)
            
            pp = props{p};
            
            switch pp
%                 case {'lh','a'}
%                     OUT.totmf.(pp) = zeros(size(bins.(pp)));
%                     OUT.ncounts.(pp) = zeros(size(bins.(pp)));
%                 otherwise
                case {'atot','atotn'}
                    OUT.(pp) = zeros(60,80,50);
                case {'adist','adistn'}
                    OUT.(pp) = zeros(size(XI));
                case {'dmedian','d25','d75','d15','d85'}
                    OUT.(pp) = zeros(80,length(files{q})); 
                case {'ad','mfd'}
                    OUT.(pp) = [];
                otherwise
                    OUT.(pp) = nan(1,length(files{q}));
            end
            
        end
        
        files{q}
        
        % NOW LOAD FILES
        for i = 1:length(files{q})
            
            disp(files{q}(i).name)
            try
                load([direcs{q} files{q}(i).name])
            catch err
                err
                continue
            end
            
            switch q
                case 1
                    Q = Airs;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % fix outliers near granule edge:
                    region = D2 > 350;
                    Q.ST.HA(region & Q.ST.HA > 4) = NaN;
%                     nanlocs = isnan(Q.ST.HA);
%                     Q.ST.HA(nanlocs) = 0;
%                     smoo = smoothn(Q.ST.HA,[11 21 5]);
%                     region = D2 > 350;
%                     Q.ST.HA(region) = smoo(region);
%                     Q.ST.HA(nanlocs) = NaN;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case {2,3}
                    Q = Model;
            end
            
%             % force WESTWARD assumption (useful for meridional flux)
%             posinds = Q.ST.F2 > 0;
%             Q.ST.F1(posinds) = -Q.ST.F1(posinds);
%             Q.ST.F3(posinds) = -Q.ST.F3(posinds);
%             Q.ST.F2          = -abs(Q.ST.F2);
%             Q.MF1(posinds)   = -Q.MF1(posinds);
%             Q.MF2            = -abs(Q.MF2);

            % try UPWARD assumption instead?
            posinds = Q.ST.F3 > 0;
            Q.ST.F1(posinds) = -Q.ST.F1(posinds);
            Q.ST.F2(posinds) = -Q.ST.F2(posinds);
            Q.ST.F3          = -abs(Q.ST.F3);
            Q.MF1(posinds)   = -Q.MF1(posinds);
            Q.MF2(posinds)   = -Q.MF2(posinds);
            
            amp   = Q.ST.HA(inds) ./ sf(inds);
            kh = quadadd(Q.ST.F1,Q.ST.F2);
            
            % extract amplitude and standard deviation:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            W.a = nanmean(amp(:));
            W.astd = nanstd(amp(:));
            
            % extract median amplitudes and percentiles:
            W.amedian = nanmedian(amp(:));
            for c = 1:length(percentiles)
                W.(['p' num2str(percentiles(c))]) = prctile(amp(:),percentiles(c));
            end
            
%             % total amplitudes for mean...
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             W.atot =  Q.ST.HA ./ sf;
%             W.atotn = double(~isnan(Q.ST.HA));
            
            % distance from island?
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for this, just stack all the amplitudes together and we'll
            % take the median and percentiles afterwards
            W.ad = Q.ST.HA ./ sf;
            W.ad = nanmean(W.ad(:,:,inrange(Airs.Alt,zlims)),3)';
            W.ad = W.ad(:,abs(d1) <= r);
            
            W.mfd = quadadd(Q.MF1,Q.MF2);
            W.mfd = nanmean(W.mfd(:,:,inrange(Airs.Alt,zlims)),3)';
            W.mfd = W.mfd(:,abs(d1) <= r);
            
            % extract directional mf
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            mf = quadadd(Q.MF1,Q.MF2);
            mf = mf(inds);
            mf1 = linearise(Q.MF1(inds));
            mf2 = linearise(Q.MF2(inds));
            
            W.mf = nanmean(mf);
            W.mf1 = nanmean(mf1);
            W.mf2 = nanmean(mf2);
            
            mf1pos = mf1;
            mf1neg = mf1;
            mf2pos = mf2;
            mf2neg = mf2;
            
            mf1pos(mf1 < 0) = 0;
            mf1neg(mf1 > 0) = 0;
            mf2pos(mf2 < 0) = 0;
            mf2neg(mf2 > 0) = 0;
            
            W.mf1pos = nanmean(mf1pos);
            W.mf1neg = nanmean(mf1neg);
            W.mf2pos = nanmean(mf2pos);
            W.mf2neg = nanmean(mf2neg);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            W.Time = datenum(Q.Time);
            
            for p = 1:length(props)
                
                pp = props{p};
                
                switch pp
%                     case {'lh','a'}
%                         
%                         % bin2mat!
%                         totmf.(pp)    = bin2mat(W.(pp),ones(size(W.(pp))),quadadd(Q.MF1(inds),Q.MF2(inds)),XI.(pp),YI.(pp),@nansum);
%                         ncounts.(pp)  = bin2mat(W.(pp),ones(size(W.(pp))),ones(size(W.(pp))),XI.(pp),YI.(pp),@sum);
%                         
%                         % and subscribe only where there was data
%                         nanlocs = isnan(totmf.(pp)) | isnan(ncounts.(pp));
%                         OUT.totmf.(pp)(~nanlocs) = OUT.totmf.(pp)(~nanlocs) + totmf.(pp)(~nanlocs);
%                         OUT.ncounts.(pp)(~nanlocs) = OUT.ncounts.(pp)(~nanlocs) + ncounts.(pp)(~nanlocs);
%                         

                    case {'atot','atotn','adist','adistn'}
                        
                        OUT.(pp) = OUT.(pp) + W.(pp);
                        
                    case {'dmedian','d25','d75','d15','d85'}
                        
                        OUT.(pp)(:,i) = W.(pp);
                        
                    case {'ad','mfd'}
                        
                        OUT.(pp) = cat(2,OUT.(pp),W.(pp));
                        
                    otherwise
                        
                        OUT.(pp)(i) = W.(pp);
                        
                end
                
                
            end
            
            
        end
        
        % finally, assign output:
        switch q
            case 1
                A = OUT;
            case 2
                M = OUT;
            case 3
                MA = OUT;
        end
        
        
    end
    
    
    
    return
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADD WIND AGAINST TIME?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading Wind...')

windfiles = [...
    dir([winddirec '201306*.mat']) ; ...
    dir([winddirec '201307*.mat']) ; ...
    dir([winddirec '201506*.mat']) ; ...
    dir([winddirec '201507*.mat']) ];

U = nan(118,length(windfiles));
V = nan(118,length(windfiles));
windtime = nan(1,length(windfiles));

for i = 1:length(windfiles)
    
    load([winddirec windfiles(i).name]);
    
    U(:,i) = double(Wind.u);
    V(:,i) = double(Wind.v);
    windtime(i) = double(Wind.Time);
    
    zz = double(Wind.z./1000);
    
end

new_zrange = 0:1:70;

[X,Y] = meshgrid(1:length(windfiles),zz);
[XI,YI] = meshgrid(1:length(windfiles),new_zrange);

% Now interpolate to something useful:
Ui = interp2(X,Y,U,XI,YI);
Vi = interp2(X,Y,V,XI,YI);

% Split into years:
dv = datevec(windtime);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % or load saved version
else
    
%     load('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_TimeSeries_MeridMF')
    load('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_TimeSeries_MeridMF_4dp')
    
end





return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_TimeSeries_MeridMF','A','M','MA','r','U','V','Ui','Vi','windtime','new_zrange','dv')
% disp('Saved.')

save('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_TimeSeries_MeridMF_4dp','A','M','MA','r','U','V','Ui','Vi','windtime','new_zrange','dv')
disp('Saved.')


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTY WOTTY PLOT PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.9 1])

%-------------------------------------------------------
vert_gap    = 0.08;    horz_gap    = 0.06;
lower_marg  = 0.06;    upper_marg  = 0.05;
left_marg   = 0.085;     right_marg  = 0.065;

rows = 3; cols = 2;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 18;

%--------------------------------------------------------

prop = 'MF2'; % A | MF1 | MF2 | KH | KZ | Var

xlims = {...
    [datenum(2013,07,01,16,00,00) datenum(2013,07,31,17,00,00)],...
    [datenum(2015,06,12,16,00,00) datenum(2015,07,07,03,00,00)]};

nclev = 120;

tricolor = [mcolor('red') ; [0.3 0.3 0.3] ; mcolor('blue')];

tbpx    = {...
    A.Time,A.Time,...
    M.Time,M.Time,...
    MA.Time,MA.Time};

% % % % trim to selected time regions?
% % % tbpx    = {...
% % %     A.Time,A.Time;...
% % %     M.Time,M.Time;...
% % %     MA.Time,MA.Time};
% % % for i = 1:2
% % %     tbpx{1,i}(~inrange(tbpx{1,i},xlims{i})) = NaN;
% % %     tbpx{2,i}(~inrange(tbpx{2,i},xlims{i})) = NaN;
% % %     tbpx{3,i}(~inrange(tbpx{3,i},xlims{i})) = NaN;
% % % end
% % % tbpx = reshape(tbpx',1,6);

years = {2013,2015,2013,2015,2013,2015};

leglabs = {...
    30:10:50,30:10:50,...
    20:10:60,20:10:60,...
    20:10:60,20:10:60};

tits = {'AIRS','AIRS',' Model',' Model',' Model as AIRS',' Model as AIRS'};

switch prop
    case 'A'
        ylims = {...
            [0 8],[0 8],...
            [0 20],[0 20],...
            [0 8],[0 8]};
        ytix = {0:2:10,0:2:10,...
            0:4:20,0:4:20,...
            0:10,0:10};
        yminortix = {0:0.5:10,0:0.5:10,...
            0:1:20,0:1:20,...
            0:0.5:10,0:0.5:10};
        ylab = 'Wave Amplitude (K)';
    case 'MF1'
        ylims = {...
            [-15 10],[-15 10],...
            [-300 100],[-300 100],...
            [-15 10],[-15 10]};
        ytix = {-50:5:50,-50:5:50,...
            -300:100:100,-300:100:100,...
            -50:5:50,-50:5:50};
        yminortix = {-15:10,-15:10,...
            -300:25:100,-300:25:100,...
            -15:10,-15:10};
        ylab = '{\bf{MF_y}} (mPa)';
    case 'MF2'
        ylims = {...
            [-20 5],[-20 5],...
            [-400 100],[-400 100],...
            [-20 5],[-20 5]};
        ytix = {-50:5:50,-50:5:50,...
            -400:100:100,-400:100:100,...
            -50:5:50,-50:5:50};
        yminortix = {-30:10,-30:10,...
            -400:25:100,-400:25:100,...
            -30:10,-30:10};
        ylab = '{\bf{MF_x}} (mPa)';
end


colors = {...
    mcolor('red') ; mcolor('red') ; ...
    [0.3 0.3 0.3] ; [0.3 0.3 0.3] ; ...
    mcolor('blue') ; mcolor('blue')};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ax = 1:6
    
    yr = years{ax};
    yrstr = num2str(yr);
    
    subplot(rows,cols,ax)
    
    switch ax
        case {1,2}
            Q = A;
        case {3,4}
            Q = M;
        case {5,6}
            Q = MA;
    end
    
    % set some timeinds?
    switch ax
        case 1 % AIRS
            timeinds = inrange(Q.Time,xlims{1} + [0 0]);
        case 2 % AIRS
            timeinds = inrange(Q.Time,xlims{2} + [0 0]);
        case {3,5} % MODELS
            timeinds = inrange(Q.Time,xlims{1} + [0 0]);
        case {4,6} % MODELS
            timeinds = inrange(Q.Time,xlims{2} + [0 0]);
    end
%     timeinds = ones(size(Q.Time)) == 1;
        
    % quick smooth?
    smoo = 5;
    
    switch ax
        case {1,2} % AIRS
            % don't smooth, it's too sparse
            switch prop
                case 'MF1'
                    mf = Q.mf1(timeinds);
                    mfpos = Q.mf1pos(timeinds);
                    mfneg = Q.mf1neg(timeinds);
                case 'MF2'
                    mf = Q.mf2(timeinds);
                    mfpos = Q.mf2pos(timeinds);
                    mfneg = Q.mf2neg(timeinds);
            end
            
        case {3,4,5,6} % models
            
            % smooth models
            switch prop
                case 'MF1'
                    mf = movmean(Q.mf1(timeinds),smoo,1,'omitnan');
                    mfpos = movmean(Q.mf1pos(timeinds),smoo,1,'omitnan');
                    mfneg = movmean(Q.mf1neg(timeinds),smoo,1,'omitnan');
                    Q.amedian = movmean(Q.amedian(timeinds),smoo,2,'omitnan');
                case 'MF2'
                    mf = movmean(Q.mf2(timeinds),smoo,1,'omitnan');
                    mfpos = movmean(Q.mf2pos(timeinds),smoo,1,'omitnan');
                    mfneg = movmean(Q.mf2neg(timeinds),smoo,1,'omitnan');
                    Q.amedian = movmean(Q.amedian(timeinds),smoo,2,'omitnan');
            end
    end
    
%     mfpos(mfpos > 0.005) = 0.005;
    
    switch prop
        case {'MF','MF1','MF2'}
            % filled MF:
            hold on; F = fill([Q.Time(timeinds) reverse(Q.Time(timeinds))],1000. * [mfpos reverse(mfneg)],'k');
            F.FaceColor = rgbtrip(0);
            F.FaceAlpha = 0.15;
            F.LineStyle = 'none';
        case 'A'
            switch ax
                case {1,2}
                    fillline1 = [Q.p5(timeinds) reverse(Q.p95(timeinds))];
                    fillline2 = [Q.p15(timeinds) reverse(Q.p85(timeinds))];
                    fillline3 = [Q.p25(timeinds) reverse(Q.p75(timeinds))];
                otherwise
                    fillline1 = [movmean(Q.p5(timeinds),smoo,2,'omitnan') reverse(movmean(Q.p95(timeinds),smoo,2,'omitnan'))];
                    fillline2 = [movmean(Q.p15(timeinds),smoo,2,'omitnan') reverse(movmean(Q.p85(timeinds),smoo,2,'omitnan'))];
                    fillline3 = [movmean(Q.p25(timeinds),smoo,2,'omitnan') reverse(movmean(Q.p75(timeinds),smoo,2,'omitnan'))];
            
            end
            % filled Amplitude percentiles:
            hold on; F = fill([Q.Time(timeinds) reverse(Q.Time(timeinds))],fillline1,'k');
            F.FaceColor = rgbtrip(0);
            F.FaceAlpha = 0.15;
            F.EdgeAlpha = 0.25;
            F.LineWidth = 0.5;
%             F.LineStyle = 'none';
            % filled Amplitude percentiles:
            hold on; F = fill([Q.Time(timeinds) reverse(Q.Time(timeinds))],fillline2,'k');
            F.FaceColor = rgbtrip(0);
            F.FaceAlpha = 0.15;
            F.EdgeAlpha = 0.25;
            F.LineWidth = 0.5;
%             F.LineStyle = 'none';
            % filled Amplitude percentiles:
            hold on; F = fill([Q.Time(timeinds) reverse(Q.Time(timeinds))],fillline3,'k');
            F.FaceColor = rgbtrip(0);
            F.FaceAlpha = 0.15;
            F.EdgeAlpha = 0.25;
            F.LineWidth = 0.5;
%             F.LineStyle = 'none';

    end
    
    % marker size:
    msize = 24;
    
    % data lines!
    switch prop
        case 'A'
            % AIRS
            switch ax
                case {1,2}
                    % dotted zero line:
                    hold on; plot(Q.Time(timeinds),zeros(size(Q.Time(timeinds))),'linewi',1,'linest','--','color',rgbtrip(.15));
                    % data itself
                    hold on; plot(Q.Time(timeinds),Q.amedian(timeinds),'color','w','linewi',3.25,'marker','.','markersize',msize+4); grid on;
                    hold on; plot(Q.Time(timeinds),Q.amedian(timeinds),'color',colors{ax},'linewi',2.5,'marker','.','markersize',msize); grid on;
                otherwise
                    % MODELS
                    % dotted zero line:
                    hold on; plot(Q.Time(timeinds),zeros(size(Q.Time(timeinds))),'linewi',1,'linest','--','color',rgbtrip(.15));
                    % data itself
                    hold on; plot(Q.Time(timeinds),Q.amedian(timeinds),'color','w','linewi',3.25); grid on;
                    hold on; plot(Q.Time(timeinds),Q.amedian(timeinds),'color',colors{ax},'linewi',2.5); grid on;
            end
            % MODELS
        otherwise
            % MF
            switch ax
                case {1,2}
                    % fill symbols
                    hold on; plot(Q.Time(timeinds),1000*mfpos,'color',rgbtrip(.8),'linest','none','marker','.','markersize',16); grid on;
                    hold on; plot(Q.Time(timeinds),1000*mfneg,'color',rgbtrip(.8),'linest','none','marker','.','markersize',16); grid on;
                    % dotted zero line:
                    hold on; plot(Q.Time(timeinds),zeros(size(Q.Time(timeinds))),'linewi',1,'linest','--','color',rgbtrip(.15));
                    % data itself
                    hold on; plot(Q.Time(timeinds),1000*mf,'color','w','linewi',3,'marker','.','markersize',msize+4); grid on;
                    hold on; plot(Q.Time(timeinds),1000*mf,'color',colors{ax},'linewi',2.5,'marker','.','markersize',msize); grid on;
                otherwise
                    % dotted zero line:
                    hold on; plot(Q.Time(timeinds),zeros(size(Q.Time(timeinds))),'linewi',1,'linest','--','color',rgbtrip(.15));
                    % data itself
                    hold on; plot(Q.Time(timeinds),1000*mf,'color','w','linewi',3); grid on;
                    hold on; plot(Q.Time(timeinds),1000*mf,'color',colors{ax},'linewi',2.5); grid on;
            end
    end
    
    xlim(xlims{2-mod(ax,2)})
    % Dates...
    axx = gca;
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
    
    % Axes formatting
    xlabel(yrstr)
    set(gca,'tickdir','out','linewi',1.5,'clipping','off')
    ylabel(ylab)
    axx.YMinorTick = 'on';

    
    % COLORBOX LABELS
    nudgefactor = [1.025 1.025 1 1 0.975 0.975];
    switch ax
        case {1,3,5}
            
            axx2 = axes('position',axx.Position,'xlim',axx.XLim,'ylim',axx.YLim,'color','none','xcolor','none','ycolor','none');
            axx2.Clipping = 'off';
            
            axx2.YLim = ylims{ax};
            
            T = text(axx2,axx2.XLim(1)-5,axx2.YLim(1) + nudgefactor(ax)*(abs(diff(axx2.YLim))/2),tits{ax},'color','w','fontsize',1.5*fs,'rotation',90,'fontweight','bold','HorizontalAlignment','center');
            T.Layer = 'front';

            % Put a coloured box under it:
            xl = [T.Extent(1) T.Extent(1)+1.75];
            yl = axx2.YLim;
            
            hold on; P = patch(axx2,xl([1 2 2 1 1]),yl([1 1 2 2 1]),mcolor('red'));
            
            P.LineStyle = 'none';
            P.EdgeColor = 'none';
            
            
            switch ax
                case 1, P.FaceColor = tricolor(1,:); % AIRS
                case 3, P.FaceColor = tricolor(2,:); % MODEL
                case 5, P.FaceColor = tricolor(3,:); % MODEL AS AIRS
            end
            
    end
    
%     switch prop
%         case 'A'
%             switch ax
%                 case 2
%                     hold on; plot(axx.XLim,axx.YLim([2 2]),'color',[0 0 0 0.15],'linewi',1);
%             end
%     end
    
    drawnow;
    nph_text([-0.02 1.05],['(' alphabet(ax) ')'],'fontsize',1.75*fs,'HorizontalAlignment','center','textborder','w');
    drawnow;
    
    setfont(fs)
    
%     return
    
end










return


%% EXPORT? ================================================================

switch prop
    case 'A'
        savename = ['~/Desktop/AirsModel_Timeseries_MedianAmp'];
    case 'MF1'
        savename = ['~/Desktop/AirsModel_Timeseries_MeridMF'];
    case 'MF2'
        savename = ['~/Desktop/AirsModel_Timeseries_ZonalMF'];
end

disp(['Exporting to ' savename '...'])

nph_saveas(gcf,savename,'png')

disp('Done.')







return











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ALSO PLOT A WIND AGAINST HEIGHT FIGURE TO FIT IN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.9 1])

%-------------------------------------------------------
vert_gap    = 0.08;    horz_gap    = 0.06;
lower_marg  = 0.06;    upper_marg  = 0.05;
left_marg   = 0.085;     right_marg  = 0.065;

rows = 3; cols = 2;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 18;

%--------------------------------------------------------

prop = 'MF1';


xlims = {...
    [datenum(2013,07,01,16,00,00) datenum(2013,07,31,17,00,00)],...
    [datenum(2015,06,12,16,00,00) datenum(2015,07,07,03,00,00)]};

timerange.yr2013 = xlims{1};
timerange.yr2015 = xlims{2};

% timerange.yr2013 = [datenum(2013,07,1.5) datenum(2013,08,01)];
% timerange.yr2015 = [datenum(2015,06,12.5) datenum(2015,07,7)];

nclev = 120;

% Choose Colors:
colorheights    = [20 30 40 50 60];
colors          = [...
    [0 0 0] ; ...
    mcolor('purple') ; ...
    mcolor('blue') ; ...
    mcolor('light_blue') ; ...
    rgbtrip(.75)];

tricolor = [mcolor('red') ; [0.3 0.3 0.3] ; mcolor('blue')];

tbpx    = {...
    A.Time,A.Time,...
    M.Time,M.Time,...
    MA.Time,MA.Time};

years = {2013,2015,2013,2015,2013,2015};

leglabs = {...
    30:10:50,30:10:50,...
    20:10:60,20:10:60,...
    20:10:60,20:10:60};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ax = 1:2
    
    yr = years{ax};
    yrstr = num2str(yr);
    
    subplot(rows,cols,ax)
    
    switch prop
        case {'MF','A'}
            U = movmean(quadadd(Ui(:,dv(:,1) == yr),Vi(:,dv(:,1) == yr)),3,2,'omitnan');
            clims = [0 130];
            cbarticks = -180:20:180;
            clev = -180:10:180;
            cmap = colormap(gca,cbrew('Blues',22));
            boxcolor = rgbtrip(.5);
            windname = 'Model Wind';
            labelcolors = {'w','k'};
        case 'MF1'
            U = movmean(Vi(:,dv(:,1) == yr),3,2,'omitnan');
            clims = [-50 50];
            cbarticks = -180:20:180;
            clev = -180:7.5:180;
            cmap = flipud(colormap(gca,cbrew('nph_BuOr',22)));
            boxcolor = nph_saturate(mcolor(2),0.75);
            windname = 'Meridional Wind';
            labelcolors = {'k','w'};
        case 'MF2'
            U = movmean(Ui(:,dv(:,1) == yr),6,2,'omitnan');
            clims = [-130 130];
            cbarticks = -200:40:200;
            clev = [-180:10:-10 -5:5:5 10:10:180];
            cmap = colormap(gca,flipud(cbrew('nph_BuOr',32)));
            boxcolor = nph_saturate(mcolor(1),0.75);
            windname = 'Zonal Wind';
            labelcolors = {'w','k'};
    end 
    
    [X,Y] = meshgrid(windtime(dv(:,1) == yr),new_zrange);
    
    % fix bottom levels
    U(isnan(U)) = 0;
    
    % quick smooth...
    tbp = imgaussfilt(U,[0.75 0.75]);
    
    % fix to axes limits:
    ylim([0 65])
    xlim(timerange.(['yr' yrstr]))
    axx = gca;
    xl = axx.XLim; yl = axx.YLim;
    tbp(~inrange(X,xl)) = NaN;
    tbp(~inrange(Y,yl)) = NaN;
    
    ylabel('Altitude (km)')
    
    %%% wind contours:
    hold on; contourf(X,Y,tbp,clev,'edgecolor','none');
    
    %%% some contour lines for merid
    switch prop
        case {'MF1','MF2'}
            hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor',[.3 .3 .3]);
    end
    colormap(gca,cmap);
    clim(clims)
    
    %     % saturate colormap?
    %     cmap_hsv = rgb2hsv(cmap);
    %     cmap_hsv(:,2) = 1.2.*cmap_hsv(:,2); cmap_hsv(cmap_hsv > 1) = 1;
    %     cmap_sat = hsv2rgb(cmap_hsv);
    %     colormap(gca,cmap_sat)
    
    % apply limits again in case contours have moved
    ylim([0 65])
    xlim(timerange.(['yr' yrstr]))
    
    
    if ax == 2
        cbar = colorbar; drawnow;
        cbar.Position = [1.032 1 1 0.6] .* cbar.Position;
        cbar.Position = [0 0.35*cbar.Position(4) 0 0] + cbar.Position;
        cbar.Ticks = cbarticks;
        cbar.Label.String = 'ms^-^1';
        cbar.TickDirection = 'out';
    end
    
    % Dates...
    axx = gca;
    axx.XTick = [...
        datenum(yr,06,[1 5:5:30]) ...
        datenum(yr,07,[  5:5:30]) ...
        datenum(yr,08,[  5:5:25])];
    axx.XMinorTick = 'on';
    axx.XAxis.MinorTickValues = datenum(yr,06,1:100);
    
    ylim([0 65])
    xlim(timerange.(['yr' yrstr]))
    datetick('x','dd mmmm','keepticks','keeplimits')
    
    % Axes formatting
    xlabel(yrstr)
    set(gca,'tickdir','out','linewi',1.5,'layer','top','clipping','off')
    
    % box around the axes
    xl = axx.XLim; yl = axx.YLim;
    hold on; plot(xl([2 2 1]),yl([1 2 2]),'color','k','linewi',1.5);
    
    % Coloured box for label
    switch ax
        case {1,3,5}
            
            T = text(axx.XLim(1)-5,0.5*axx.YLim(2),windname,'color','w','fontsize',1.5*fs,'rotation',90,'fontweight','bold','HorizontalAlignment','center');
            T.Layer = 'front';
            
            % Put a coloured box under it:
            xl = [T.Extent(1) T.Extent(1)+1.75];
            yl = get(gca,'ylim');
            
            hold on; P = patch(xl([1 2 2 1 1]),yl([1 1 2 2 1]),boxcolor);
            
            P.LineStyle = 'none';
            P.EdgeColor = 'none';
            
    end
    
    drawnow;
    nph_text([-0.02 1.05],['(' alphabet(ax+2) ')'],'fontsize',1.75*fs,'HorizontalAlignment','center','textborder','w');
    drawnow;
    
    setfont(fs)
    
    drawnow;
    
end



% return
% 
% 
% %% EXPORT? ================================================================


switch prop
    case {'MF','A'}
        savename = ['~/Desktop/AirsModel_Timeseries_absWind'];
    case 'MF1'
        savename = ['~/Desktop/AirsModel_Timeseries_MeridWind2'];
    case 'MF2'
        savename = ['~/Desktop/AirsModel_Timeseries_ZonalWind'];
end

disp(['Exporting to ' savename '...'])

nph_saveas(gcf,savename,'png')

disp('Done.')




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AMPLITUDE AND GINI AGAINST DISTANCE FROM ISLAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([1 0.45])

%-------------------------------------------------------
vert_gap    = 0.08;    horz_gap    = 0.06;
lower_marg  = 0.025;    upper_marg  = 0.05;
left_marg   = 0.05;     right_marg  = 0.05;

rows = 1; cols = 4;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 18;

%--------------------------------------------------------

prop = 'A'; % MF | A

tits = {'AIRS','AIRS',' Model',' Model',' Model as AIRS',' Model as AIRS'};

colors = {mcolor('red'),[0.3 0.3 0.3],mcolor('blue')};

switch prop
    case 'A'
%         ylims = {[0 3],[0 7],[0 3],[0.25 0.55]};
%         ytix = {0:0.5:10,0:10,0:0.5:10,0:0.05:1};
        ylims = {[0 3.5],[0 8],[0 3.5],[0.3 1]};
        ytix = {0:0.5:10,0:1:10,0:0.5:10,0:0.1:1};
        ylabs = {'Wave Amplitude (K)','Wave Amplitude (K)','Wave Amplitude (K)','Gini Coefficient'};
    case 'MF'
        ylims = {[0 25],[0 450],[0 25],[0.3 1]};
        ytix = {0:5:30,0:50:450,0:5:50,0:0.1:1};
        ylabs = {'GWMF (mPa)','GWMF (mPa)','GWMF (mPa)','Gini Coefficient'};
end


% xlims = {...
%     [-500 500],...
%     [-500 500],...
%     [-500 500],...
%     [-500 500]};
xlims = [-500 500];
% xtix = {[-800:200:800],[-800:200:800],[-800:200:800],[-800:200:800]};
xtix = {...
    [-1000:250:1000],...
    [-1000:250:1000],...
    [-1000:250:1000],...
    [-1000:200:1000]};


% for doing the percentage of MF in the final Gini axes...
frac = 0.9;
percs = 1:-0.1:0;
ginis = (((-1/frac) .* percs)+1);
    

d2 = linspace(-600,600,80);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ax = 1:4
    
    subplot(rows,cols,ax)
    
    switch ax
        case 1
            Q = A;
        case 2
            Q = M;
        case 3
            Q = MA;
            
    end
    
    switch ax
        case {1,2,3}
            % filled Amplitude percentiles:
            for pc = [25 15 5]
                
                switch prop
                    case 'A'
                        fillline = [prctile(Q.ad,pc,2) ; reverse(prctile(Q.ad,100-pc,2))];
                    case 'MF'
                        fillline = 1000*[prctile(Q.mfd,pc,2) ; reverse(prctile(Q.mfd,100-pc,2))];
                end
                
                hold on; F = fill([d2 reverse(d2)],fillline,'k');
                F.FaceColor = rgbtrip(0);
                F.FaceAlpha = 0.15;
                F.EdgeAlpha = 0.25;
                F.LineWidth = 1;
            end
            
            %%%% data lines!
            switch prop
                case 'A'
                    hold on; plot(d2,nanmedian(Q.ad,2),'color','w','linewi',3.5); grid on;
                    hold on; plot(d2,nanmedian(Q.ad,2),'color',colors{ax},'linewi',3); grid on;
                case 'MF'
                    hold on; plot(d2,1000*nanmedian(Q.mfd,2),'color','w','linewi',3.5); grid on;
                    hold on; plot(d2,1000*nanmedian(Q.mfd,2),'color',colors{ax},'linewi',3); grid on;
            end
            % logscale('y')
            
        otherwise
            % GINI COEFFICIENT OF AMP OR MF!
%             switch prop
%                 case 'A'
%                     hold on; plot(d2,ginicoeff(M.ad,2),'color',colors{2},'linewi',3); grid on;
%                     hold on; plot(d2,ginicoeff(A.ad,2),'color',colors{1},'linewi',3); grid on;
%                     hold on; plot(d2,ginicoeff(MA.ad,2),'color',colors{3},'linewi',3); grid on;
%                 case 'MF'
                    hold on; plot(d2,ginicoeff(M.mfd,2),'color',colors{2},'linewi',3); grid on;
                    hold on; plot(d2,ginicoeff(A.mfd,2),'color',colors{1},'linewi',3); grid on;
                    hold on; plot(d2,ginicoeff(MA.mfd,2),'color',colors{3},'linewi',3); grid on;
%             end
            
    end
    
    axx = gca;
    
%     xlim(xlims{ax})
    xlim(xlims)
    axx.XTick = xtix{ax};
    axx.XMinorTick = 'on';
    axx.XAxis.MinorTickValues = -800:50:800;
    
    ylim(ylims{ax})
    axx.YTick = ytix{ax};
    axx.YMinorTick = 'on';
    
    switch ax
        case {1,3}
            axx.YAxis.MinorTickValues = 0:0.25:10;
        case 2
            axx.YAxis.MinorTickValues = 0:0.5:10;
        case 4
            axx.YMinorGrid = 'on';
            axx.YAxis.MinorTickValues = ginis;
    end
    
    yl = ylabel(ylabs{ax});
    yl.Position = [1.15 1 1] .* yl.Position;
    
    xlabel('x (km)')
    
    axis square
    
    switch ax
        case 4
%             switch prop
%                 case 'A'
%                     ti = title('Intermittency','fontsize',1.25*fs);
%                     ti.Position(2) = 1.04*axx.YLim(2);
%                 case 'MF'
                    ti = title({'Intermittency','of GWMF'},'fontsize',1.25*fs);
                    ti.Position(2) = 1.02*axx.YLim(2);
%             end

    % Extra axes for percentage of events containing 90% of total MF
    yyaxis('right');
    axx = gca;
    axx.YLim = ylims{ax};
    axx.YColor = 'k';
    axx.YTick = ginis;
%     axx.YTickLabel = cellstr([num2str(100.*percs') repmat('%',1,length(percs))']);
    axx.YTickLabel = cellstr([num2str(100.*percs')]);
    axx.YLabel.String = {'Percentage of events','containing 90% of total GWMF'};
    
%     axx.YGrid = 'on';
%     axx.GridLineStyle = ':';
    
    end
    
    % PLOT TOPOGRAPHY
    xd = xdist([1:end 1]);
    tp = topo([1:end 1])./1000;
    % normalise to axes
    tp = tp ./ max(tp);
    tp = (tp .* 0.04 * abs(diff(axx.YLim))) + axx.YLim(1);
    hold on; fill(xd,tp,'k','edgecolor',[.9 .9 .9],'linewidth',1.5);
%     hold on; plot(xlims{ax},axx.YLim([1 1]),'color','k','linewi',1.5)
    hold on; plot(xlims,axx.YLim([1 1]),'color','k','linewi',1.5)
    
    % Axes formatting
    set(gca,'tickdir','out','linewi',1.5,'clipping','on')
    
    drawnow;
    nph_text([0.015 0.865],['(' alphabet(ax) ')'],'fontsize',1.75*fs,'HorizontalAlignment','center','textborder','w');
    drawnow;
    
    setfont(fs)
    
%     return
    
end










return


%% EXPORT? ================================================================

% also known as amp_gini_distance.png

switch prop
    case 'A'
        savename = ['~/Desktop/AirsModel_AmpGiniDistance_' prop];
    case 'MF'
        savename = ['~/Desktop/AirsModel_AmpGiniDistance_' prop];        
end

disp(['Exporting to ' savename '...'])

nph_saveas(gcf,savename,'png')

disp('Done.')







return






























