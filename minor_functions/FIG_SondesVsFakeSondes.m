
% NPH 14-Feb-2018 - revised late 2019!!
%
% Plot vertical differences and standard deviation of the difference
% between the sondes and model (fake sondes).
%
%
%

%% DIRECTORIES ============================================================

matlabdirec = '/Users/neil/Drive/MATLAB/';
sondedirec = '/Volumes/SDBlue/data/sgwex/sondes/';
fakesondedirec = '/Volumes/SDBlue/data/sgwex/fakesondes/';
winddirec = '/Volumes/SDBlue/data/sgwex/VertWindProfiles/';

% savedirec = '/Users/neil/data/sgwex/fakesondes/';

%% SELECT CAMPAIGN ========================================================

camp = 'winter';

switch camp
    case {'winter' 'w'}
        fakesondefiles02s = [...
            dir([fakesondedirec '201506*02s*FakeSonde*.mat']) ; ...
            dir([fakesondedirec '201507*02s*FakeSonde*.mat']) ];
        fakesondefiles10s = [...
            dir([fakesondedirec '201506*10s*FakeSonde*.mat']) ; ...
            dir([fakesondedirec '201507*10s*FakeSonde*.mat']) ];
    case {'summer' 's'}
        fakesondefiles02s = dir([fakesondedirec '201501*02s*FakeSonde*.mat']);
        fakesondefiles10s = dir([fakesondedirec '201501*10s*FakeSonde*.mat']);
end

%% FIND UNIQUE LAUNCHES ===================================================
% sometimes there is only a 2s or 10s launch, due to the software crashing
% or the file being corrupt, so find all the unique launches.

% we want all the 2s launches, but not all the 10s ones:
a02s = split({fakesondefiles02s.name},'_');
a02s = a02s(:,:,1);

a10s = split({fakesondefiles10s.name},'_');
a10s = a10s(:,:,1);

% check if there were any 10s launches for which we DO NOT have a 2s
% launch:
extra_10s_launches = zeros(size(a10s));
for i = 1:length(a10s)
    if ~any(strcmpi(a10s{i},a02s))
        %         disp(['got unique ' a10s{i}])
        extra_10s_launches(i) = 1;
    end
end

% select only the extra launches and add them to our file structure:
fakesondefiles = [fakesondefiles02s ; fakesondefiles10s(logical(extra_10s_launches))];

tensectwosec = [repmat(2,length(fakesondefiles02s),1)  ; repmat(10,length(fakesondefiles02s),1)];

% now put them in chronological order:
a = split({fakesondefiles.name},'_');
a = a(:,:,1);
a = datenum(a,'yyyymmdd-HHMM');
[~,chronorder] = sort(a);

fakesondefiles = fakesondefiles(chronorder);
tensectwosec = tensectwosec(chronorder);


%% LOAD SONDES AND EXTRACT VARIABLES ==================================

disp('Loading Sondes...')

desired_alt_range = linspace(0,40,400); % km

Su = nan(length(fakesondefiles),length(desired_alt_range));
Sv = nan(length(fakesondefiles),length(desired_alt_range));
St = nan(length(fakesondefiles),length(desired_alt_range));

FSu = nan(length(fakesondefiles),length(desired_alt_range));
FSv = nan(length(fakesondefiles),length(desired_alt_range));
FSt = nan(length(fakesondefiles),length(desired_alt_range));

LaunchTimes = nan(1,length(fakesondefiles));
MaxAlts =   nan(1,length(fakesondefiles));

for i = 1:length(fakesondefiles)
    
    %%%% ignore some fakesondes, no model data:
    if strcmpi(fakesondefiles(i).name,'20150109-1105_10s_FakeSonde.mat')
        continue
    end
    
    % loady load...
    load([fakesondedirec fakesondefiles(i).name])
    
    load([sondedirec fakesondefiles(i).name(1:18) 'Sonde.mat'])
    
    if any(regexp(fakesondefiles(i).name,'10s'))
        Sonde.WindSpeed = Sonde.WindSpeed .* 0.51; % convert knots to m/s
    end
    
    % project sonde wind into u and v:
    Sonde.u = Sonde.WindSpeed.*sind(Sonde.WindDirection);
    Sonde.v = Sonde.WindSpeed.*cosd(Sonde.WindDirection);
    
    % find all the non-nan, unique indices
    [a,b] = unique(FakeSonde.Alt);
    fs_des_inds = b(~isnan(a));
    
    % interpolate onto the same vertical grid:
    FakeSonde.ui = interp1(FakeSonde.Alt(fs_des_inds),FakeSonde.u(fs_des_inds),desired_alt_range,'linear');
    FakeSonde.vi = interp1(FakeSonde.Alt(fs_des_inds),FakeSonde.v(fs_des_inds),desired_alt_range,'linear');
    FakeSonde.ti = interp1(FakeSonde.Alt(fs_des_inds),FakeSonde.Temp(fs_des_inds),desired_alt_range,'linear');
    
    FSu(i,:) = FakeSonde.ui;
    FSv(i,:) = FakeSonde.vi;
    FSt(i,:) = FakeSonde.ti;
    
    % find all the non-nan, unique indices
    [a,b] = unique(Sonde.Alt);
    s_des_inds = b(~isnan(a));
    
    Sonde.ui = interp1(Sonde.Alt(s_des_inds),Sonde.u(s_des_inds),desired_alt_range,'linear');
    Sonde.vi = interp1(Sonde.Alt(s_des_inds),Sonde.v(s_des_inds),desired_alt_range,'linear');
    Sonde.ti = interp1(Sonde.Alt(s_des_inds),Sonde.Temp(s_des_inds),desired_alt_range,'linear');
    
    Su(i,:) = Sonde.ui;
    Sv(i,:) = Sonde.vi;
    St(i,:) = Sonde.ti;
    
    LaunchTimes(i) = datenum(Sonde.StartTime);
    MaxAlts(i) = nanmax(Sonde.Alt);
    
end

% try fixing wind artefacts by using missing temperature data?
nanlocs = isnan(St) | isnan(FSt);
Su(nanlocs) = NaN;
% FSv(nanlocs) = NaN;
Sv(nanlocs) = NaN;
% FSv(nanlocs) = NaN;
St(nanlocs) = NaN;

% nice, that seems to have improved things. Now try using a max difference
% argument:
du = [zeros(size(Su,1),1) abs(diff(Su,[],2))];
dv = [zeros(size(Su,1),1) abs(diff(Sv,[],2))];
difflim = 17.5;
badlocs = du > difflim | dv > difflim;

Su(badlocs) = NaN;
% FSv(badlocs) = NaN;
Sv(badlocs) = NaN;
% FSv(badlocs) = NaN;
St(badlocs) = NaN;

% FIX SPECIFIC PROBLEMS BY HAND:

Su(1,inrange(desired_alt_range,27,30)) = NaN;
Su(15,inrange(desired_alt_range,0,1)) = NaN;
Su(29,inrange(desired_alt_range,26,30)) = NaN;
Su(29,inrange(desired_alt_range,26,30)) = NaN;
Su(33,inrange(desired_alt_range,4,12)) = NaN;
Su(42,inrange(desired_alt_range,2,8)) = NaN;
Su(49,inrange(desired_alt_range,23,30)) = NaN;
Su(54,inrange(desired_alt_range,3,4)) = NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Also load some wind:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading Wind...')
windfiles = [dir([winddirec '201506*.mat']) ; dir([winddirec '201507*.mat'])];

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

[X,Y] = meshgrid(windtime,zz);
[XI,YI] = meshgrid(windtime,0:1:60);

% Now interpolate to something useful:
Ui = interp2(X,Y,U,XI,YI);
Vi = interp2(X,Y,V,XI,YI);

UU = quadadd(Ui,Vi);


% FIND STANDARD DEVIATIONS ================================================

% smooth sondes slightly before finding difference?
Su = movmean(Su,5,2,'omitnan');
Sv = movmean(Sv,5,2,'omitnan');

% get differences, SONDE - MODEL
vdiff = Sv - FSv;
udiff = Su - FSu;
tempdiff = St - FSt;

% exclude some sondes from the STD:
%     vdiff(14:19,:) = NaN; % they're full of spikes which affect the STD

vstd = std(vdiff,'omitnan');
ustd = std(udiff,'omitnan');
tempstd = std(tempdiff,'omitnan');

%get means
vmean = nanmean(vdiff);
umean = nanmean(udiff);
tempmean = nanmean(tempdiff);

% get rid of nans:
notnans = intersect(intersect(find(~isnan(vstd)),find(~isnan(ustd))),find(~isnan(tempstd)));
z = desired_alt_range(notnans);

vstd = vstd(notnans);
ustd = ustd(notnans);
tempstd = tempstd(notnans);

vmean = vmean(notnans);
umean = umean(notnans);
tempmean = tempmean(notnans);

% Smooth the std?
window = 11;
vstd = movmean(vstd,window);
ustd = movmean(ustd,window);
tempstd = movmean(tempstd,window);



%% PLOTTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.6 1])

%-------------------------------------------------------
vert_gap = 0.025;       horz_gap = 0.085;
lower_marg = 0.075;   upper_marg = 0.05;
left_marg = 0.1;   right_marg = 0.08;

rows = 9; cols = 3;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 17;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% NOW SONDES/FAKE SONDES WINNDS AND DIFFERENCES

TBP = {...
    Su,FSu,udiff,...
    Sv,FSv,vdiff};

% axpositions = {[4 7],[5 8],[6 9],[10 13],[11 14],[12 15]};
% axpositions = {[7 10],[8 11],[9 12],[13 16],[14 17],[15 18]};
axpositions = {[10 13 16],[11 14 17],[12 15 18],[19 22 25],[20 23 26],[21 24 27]};

% linespec = {'color',[.6 .6 .6],'linewi',1};

linecolors = {...
    mcolor(1),mcolor(1),mcolor(1),...
    mcolor(2),mcolor(2),rgbtrip(.6)};

xlims = {...
    [-60 100],[-60 100],[-40 40],...
    [-50 50],[-50 50],[-40 40]};


for ax = 1:6
    
    subplot(rows,cols,axpositions{ax});
    
    % plot all sondes
    if ax ~= 3 && ax ~= 6
        hold on; plot(TBP{ax}',desired_alt_range,'color',linecolors{ax},'linewi',1.25); grid on
    end
    
    % plot mean and std
    if ax == 3 || ax == 6
        
        % trim top
        if ax == 3
            TBP{ax}(:,desired_alt_range > 31.85) = NaN;
        end
        
        % STD FILL
        stdvec = nanstd(TBP{ax},1);
        meanfill = nanmean(TBP{ax},1);
        meanfill(isnan(meanfill)) = 0;
        stdvec(isnan(stdvec)) = 0;
        
        % 2 std
        stdalt = [desired_alt_range fliplr(desired_alt_range)];
        stdfill = [meanfill+(2*stdvec) fliplr(meanfill-(2*stdvec))];
        hold on; P2 = fill(stdfill,stdalt,'k');
        P2.FaceColor = rgbtrip(.8);
        P2.FaceAlpha = 1;
        P2.LineStyle = 'none';
        
%         % zero line
%         hold on; plot([0 0],[0 100],'color',rgbtrip(0),'linewi',2,'linest','-');

        % 1 std
        stdalt = [desired_alt_range fliplr(desired_alt_range)];
        stdfill = [meanfill+stdvec fliplr(meanfill-stdvec)];
        hold on; P1 = fill(stdfill,stdalt,'k');
        P1.FaceColor = rgbtrip(.5);
        P1.FaceAlpha = 1;
        P1.LineStyle = 'none';
        
        % another zero line
        hold on; plot([0 0],[0 100],'color',rgbtrip(0),'linewi',1,'linest','--');
        
        % MEAN
        hold on; plot(nanmean(TBP{ax},1),desired_alt_range,'color',rgbtrip(.1),'linewi',3); grid on

    end
    
    xlim(xlims{ax})
    ylim([0 30])
    
    axx = gca;
    
    % fill in edge of axes grid line:
    xl = axx.XLim; yl = axx.YLim;
    hold on; plot(xl([2 2]),yl([1 2]),'linewi',1.5,'color',[0 0 0 0.15]);
    
    % ticks
    axx.XTick = -100:20:100;
    axx.XMinorTick = 'on';
    axx.XAxis.MinorTickValues = -100:10:100;
    
    if ax < 3
        axx.XTick = -80:40:80;
        axx.XAxis.MinorTickValues = -100:20:100;
    end
    
    if ax <= 3
        axx.Position =  [1 1.075 1 1] .* axx.Position;
    else
        xlabel('Wind Speed (ms^-^1)')
    end
    
%     if ax == 1 || ax == 4
        ylabel('Altitude (km)')
%     end
    
    switch ax
        case 1, nph_text([0.45 1.05],upper('Sondes'),'fontsize',1.2*fs,'fontweight','bold','horizontalalignment','center');
        case 2, nph_text([0.45 1.05],upper('Model-as-Sondes'),'fontsize',1.2*fs,'fontweight','bold','horizontalalignment','center');
        case 3, nph_text([0.45 1.05],upper('Difference'),'fontsize',1.2*fs,'fontweight','bold','horizontalalignment','center');
    end
    
    switch ax
        case 1
            text(-116,17.5,'ZONAL WIND','fontsize',1.25*fs,'fontweight','bold','rotation',90,'HorizontalAlignment','center');
        case 4
            text(-85,17.5,'MERIDIONAL WIND','fontsize',1.25*fs,'fontweight','bold','rotation',90,'HorizontalAlignment','center');
    
    end
    
    drawnow;
    
    nph_text([0.03 0.875],['(' alphabet(ax+1) ')'],'fontsize',1.5*fs,'facealpha',0,'textborder','w');
    
    set(gca,'linewi',1.5,'tickdir','out')
    
    setfont(fs)
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% STEM PLOT SHOWING LAUNCH TIMES AND ALTITUDES

subplot(rows,cols,[1 2 3 4 5 6])

%%% wind speed contours:
hold on; contourf(XI,YI,UU,-20:10:180,'edgecolor','none');
colormap(gca,cbrew('Blues',21))
clim([0 100])
cbar = colorbar; drawnow;
cbar.Position = [1.04 1 1 0.6] .* cbar.Position;
cbar.Position = [0 0.35*cbar.Position(4) 0 0] + cbar.Position;
cbar.Ticks = 0:20:100;
cbar.Label.String = 'Wind Speed (ms^-^1)';

stemspec = {'color','k','linest','-','linewi',2,...
    'markerfacecolor','k','markeredgecolor','k','markersize',8};

hold on; stem(LaunchTimes,MaxAlts,stemspec{:},'color','w','markeredgecolor','w','markersize',9,'linewi',2.5);
hold on; stem(LaunchTimes,MaxAlts,stemspec{:});

a1 = gca;
axx = gca;
axx.XGrid = 'off';
axx.Layer = 'top';

axx.XTick = [datenum(2015,06,0:5:25) datenum(2015,07,0:5:30)];
axx.XAxis.MinorTickValues = datenum(2015,06,0:60);
axx.XMinorTick = 'on';
datetick('x','dd mmmm','keepticks','keeplimits')

xlim([datenum(2015,06,13) datenum(2015,07,6.6)])
ylim([0 40])

ylabel('Altitude (km)')

nph_text([-0.03 0.86],'(a)','fontsize',1.5*fs,'color','w','facealpha',0,'fitboxtotext','on','textborder','k');

set(gca,'linewi',1.5,'tickdir','out')

setfont(fs)



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT TEMPERATURE INSTEAD?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.6 1])

%-------------------------------------------------------
vert_gap = 0.025;       horz_gap = 0.085;
lower_marg = 0.075;   upper_marg = 0.05;
left_marg = 0.1;   right_marg = 0.08;

rows = 9; cols = 3;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 17;

%--------------------------------------------------------


TBP = {...
    St,FSt,tempdiff,...
    Sv,FSv,vdiff};

% axpositions = {[4 7],[5 8],[6 9],[10 13],[11 14],[12 15]};
% axpositions = {[7 10],[8 11],[9 12],[13 16],[14 17],[15 18]};
axpositions = {[10 13 16],[11 14 17],[12 15 18],[19 22 25],[20 23 26],[21 24 27]};

% linespec = {'color',[.6 .6 .6],'linewi',1};

linecolors = {...
    mcolor(4),mcolor(4),mcolor(4),...
    mcolor(2),mcolor(2),rgbtrip(.6)};

xlims = {...
    [175 300],[175 300],[-15 15],...
    [-50 50],[-50 50],[-40 40]};

for ax = 1:3
    
    subplot(rows,cols,axpositions{ax});
    
    % plot all sondes
    if ax ~= 3 && ax ~= 6
        hold on; plot(TBP{ax}',desired_alt_range,'color',linecolors{ax},'linewi',1.25); grid on
    end
    
    % plot mean and std
    if ax == 3 || ax == 6
        
        % trim top
        if ax == 3
            TBP{ax}(:,desired_alt_range > 31.85) = NaN;
        end
        
        % STD FILL
        stdvec = nanstd(TBP{ax},1);
        meanfill = nanmean(TBP{ax},1);
        meanfill(isnan(meanfill)) = 0;
        stdvec(isnan(stdvec)) = 0;
        
        % 2 std
        stdalt = [desired_alt_range fliplr(desired_alt_range)];
        stdfill = [meanfill+(2*stdvec) fliplr(meanfill-(2*stdvec))];
        hold on; P2 = fill(stdfill,stdalt,'k');
        P2.FaceColor = rgbtrip(.8);
        P2.FaceAlpha = 1;
        P2.LineStyle = 'none';
        
%         % zero line
%         hold on; plot([0 0],[0 100],'color',rgbtrip(0),'linewi',2,'linest','-');

        % 1 std
        stdalt = [desired_alt_range fliplr(desired_alt_range)];
        stdfill = [meanfill+stdvec fliplr(meanfill-stdvec)];
        hold on; P1 = fill(stdfill,stdalt,'k');
        P1.FaceColor = rgbtrip(.5);
        P1.FaceAlpha = 1;
        P1.LineStyle = 'none';
        
        % another zero line
        hold on; plot([0 0],[0 100],'color',rgbtrip(0),'linewi',1,'linest','--');
        
        % MEAN
        hold on; plot(nanmean(TBP{ax},1),desired_alt_range,'color',rgbtrip(.1),'linewi',3); grid on

    end
    
    xlim(xlims{ax})
    ylim([0 30])
    
    axx = gca;
    
    % fill in edge of axes grid line:
    xl = axx.XLim; yl = axx.YLim;
    hold on; plot(xl([2 2]),yl([1 2]),'linewi',1.5,'color',[0 0 0 0.15]);
    
    % ticks
    switch ax
        case 3
            axx.XTick = -20:10:20;
            axx.XMinorTick = 'on';
            axx.XAxis.MinorTickValues = -20:5:20;
        otherwise
            axx.XTick = 0:25:400;
            %             axx.XMinorTick = 'on';
            %             axx.XAxis.MinorTickValues = 0:10:400;
    end
    
    %     axx.XTick = -100:20:100;
    %     axx.XMinorTick = 'on';
    %     axx.XAxis.MinorTickValues = -100:10:100;
    
%     if ax < 3
%         axx.XTick = -80:40:80;
%         axx.XAxis.MinorTickValues = -100:20:100;
%     end
    
%     if ax <= 3
%         axx.Position =  [1 1.075 1 1] .* axx.Position;
%     else
        xlabel('(K)')
%     end
    
%     if ax == 1 || ax == 4
        ylabel('Altitude (km)')
%     end
    
    switch ax
        case 1, nph_text([0.45 1.05],upper('Sondes'),'fontsize',1.2*fs,'fontweight','bold','horizontalalignment','center');
        case 2, nph_text([0.45 1.05],upper('Model-as-Sondes'),'fontsize',1.2*fs,'fontweight','bold','horizontalalignment','center');
        case 3, nph_text([0.45 1.05],upper('Difference'),'fontsize',1.2*fs,'fontweight','bold','horizontalalignment','center');
    end
    
    switch ax
        case 1
            text(135,15.5,'TEMPERATURE','fontsize',1.25*fs,'fontweight','bold','rotation',90,'HorizontalAlignment','center');
%         case 4
%             text(-85,17.5,'MERIDIONAL WIND','fontsize',1.25*fs,'fontweight','bold','rotation',90,'HorizontalAlignment','center');
    end
    
    drawnow;
    
    nph_text([0.03 0.875],['(' alphabet(ax) ')'],'fontsize',1.5*fs,'facealpha',0,'textborder','w');
    
    set(gca,'linewi',1.5,'tickdir','out')
    
    setfont(fs)
    
end











%% EXPORT? ================================================================


savename = '~/Desktop/SondesVsFakeSondes';
% savename = '~/Desktop/SondesVsFakeSondes_Temperature';

disp(['Exporting to "' savename '"...'])
nph_saveas(gcf,savename,'png')

disp('Done.')



return












