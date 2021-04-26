


% Average wave amplitudes against height for AIRS, Model and Model as AIRS.
% Show depatures from exponential growth with altitude.



runagain = 1;


%% DIRECTORIES ============================================================

matlabdirec         = '/Users/neil/Drive/MATLAB/';
airsdirec           = '/Volumes/SDBlue/data/sgwex/3DST/Airs/AG_sm/';
% modeldirec          = '/Volumes/SDBlue/data/sgwex/3DST/Model/AG/';
modelasairsdirec    = '/Volumes/SDBlue/data/sgwex/3DST/Model_as_AIRS/AG_sm/';


%% REGIONS ================================================================

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

%%%% IGNORING ALTITUDES FOR THIS ONE, APPLY IT LATER
zlims = [25 45];

% keep awayyy from the edges
edgebox = inrange(D1,pm(h/2)) & inrange(D2,pm(L/2));

% Region A - Upwind of the island
Regions.A = ~ring & edgebox & inrange(D2,[-(L/2) (L/2)-d2]);

% Region B - Over/Downwind of the island
Regions.B = (ring | inrange(D2,[(L/2)-d2 (L/2)])) & edgebox;

% Region C - Both!
Regions.C = Regions.A | Regions.B;


%% VERT PROFILE =============================================================

v3 = 1.5:1.5:75;

[D1,D2,D3] = ndgrid(...
    linspace(-450,450-15,90),...
    linspace(-600,600-15,120),...
    v3);

R = quadadd(D1,D2);

% define vert height levels
zlev = 1.5:1.5:75;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

airsfiles = dir([airsdirec '/20*.mat']);
modelasairsfiles = dir([modelasairsdirec '/20*.mat']);

amplim = 0;

props = {'airs','modelasairs','model'};

if runagain
    
    for p = 1
        
        prop = props{p};
        
        switch prop
            case 'airs'
                files = airsfiles;
                direc = airsdirec;
            case 'modelasairs'
                files = modelasairsfiles;
                direc = modelasairsdirec;
            case 'model'
                load('/Users/neil/data/upwp_vpwp_tp_average.mat')
                R = {'A','B','C'};
                H = 7;
                dens = 1.225 .* exp(-M.z/H);
                Dens = repmat(dens',1,length(M.Time));
                for r = 1:3
                    M.(R{r}).A     = M.(R{r}).tp;
                    M.(R{r}).MF    = quadadd(Dens.*M.(R{r}).upwp,Dens.*M.(R{r}).vpwp);
                end
                continue
        end
        
        sz = [90 120 50];
        
        %%%% REGIONS
        R = {'A','B','C'};
        for r = 1:3
            P.(R{r}).A    = zeros(50,length(files));
            P.(R{r}).MF    = zeros(50,length(files));
        end
        
        %%%% LOADYLOAD
        for i = 1:length(files)
            
            disp(files(i).name)
            
            try
                load([direc files(i).name])
            catch err
                err
                continue
            end
            
            switch prop
                case 'airs'
                    Q = Airs;
                case 'modelasairs'
                    Q = Model;
            end
            
%             % deal with some artefacts from the scan track edge:
%             if strcmpi(files(i).name(1:13),{'20150615-0300'})
%                 rng = {1:45,100:120,1:50};
%                 Q.Tpg(rng{:}) = NaN;
%                 Q.ST.HA(rng{:}) = NaN;
%                 Q.MF1(rng{:}) = NaN;
%                 Q.MF2(rng{:}) = NaN;
%                 %                 hold on; pcolor(Q.Tpg(:,:,27)); shat; clim([-5 5])
%                 %                 title(files(i).name(1:13))
%                 %                 drawnow;
%                 %                 pause
%             end
            
            %%%% REGIONS
            R = {'A','B','C'};
            for r = 1:3
                
                amp = Q.ST.HA;
                mf = quadadd(Q.MF1,Q.MF2);
                
                % trim to region:
                amp(~Regions.(R{r})) = NaN;
                mf(~Regions.(R{r}))  = NaN;
                
                P.(R{r}).A(:,i)  = squeeze(nanmean(amp,[1 2]));
                P.(R{r}).MF(:,i) = squeeze(nanmean(mf,[1 2]));
                
            end
            
        end
        
        
        %%%% assign
        switch prop
            case 'airs'
                A = P;
            case 'modelasairs'
                MA = P;
        end
        
    end % next dataset
    
    
else
    
    disp('Loading saved fields...')
    load('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_ExpAmp_AG.mat')
    
    return
    
    %% %%%% SAVE?
    
    save('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_ExpAmp_AG.mat','A','MA','M');
    disp('Saved.')
    
    
end



return







return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PLOTTY WOTTY PLOT PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.85 0.45])

%-------------------------------------------------------
vert_gap    = 0.085;    horz_gap    = 0.055;
lower_marg  = 0.15;    upper_marg  = 0.1;
left_marg   = 0.055;     right_marg  = 0.05;

rows = 1; cols = 3;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 18;

%--------------------------------------------------------

prop = 'A'; % A | MF | KH | KZ | Var

R = 'B';


Acol = [200 10 0]./255;
Mcol = [.2 .2 .2];
MAcol = mcolor('blue');

lw = 2.5;

zz = [0 80];
inds = inrange(zlev,[20 61.5]);


% AIRS vertical window, as applied in 3DST:
% lower_height_index = 12; % corresponds to 12*1.5 = 18km
lower_height_index = 13; % corresponds to 13*1.5 = 19.5km
upper_height_index = 40; % corresponds to 40*1.5 = 60km
tapering_length = 6; % corresponds to 6*1.5 = 9km window at edge
OUT = sgwex_airs_apply_height_window(ones(60,80,50),lower_height_index,upper_height_index,tapering_length);
winvec = squeeze(OUT(1,1,:));
window = 1 - squeeze(OUT(30,:,:));
% window(xinds) = NaN;

% make plotting vectors
smoo = 3;
A.(R).Avec   = nanmean(A.(R).A,2);
MA.(R).Avec  = movmean(nanmean(MA.(R).A,2),smoo,1,'omitnan');
M.(R).Avec   = movmean(nanmean(M.(R).A,2),smoo,1,'omitnan');

A.(R).MFvec  = 1000*nanmean(A.(R).MF,2);
MA.(R).MFvec = 1000*movmean(nanmean(MA.(R).MF,2),smoo,1,'omitnan');
M.(R).MFvec  = 1000*movmean(nanmean(M.(R).MF,2),smoo,1,'omitnan');

% set airs visible limits
A.(R).Avec(~inds)  = NaN;
A.(R).MFvec(~inds) = NaN;
MA.(R).Avec(~inds)  = NaN;
MA.(R).MFvec(~inds) = NaN;

% trim off the tropopause:
M.(R).Avec(zlev <= 10)  = NaN;
M.(R).MFvec(zlev <= 10) = NaN;

% fix model reflection for high altitude
% M.(R).Avec(zlev > 67) = 0.98*M.(R).Avec(zlev > 67);
smooo = movmean(M.(R).Avec,5,1,'omitnan');
M.(R).Avec(zlev > 64) = smooo(zlev > 64);
M.(R).Avec(zlev > 69)  = NaN;
M.(R).MFvec(zlev > 69) = NaN;

for ax = 1:2
    
    subplot(rows,cols,ax)
    
    % AXES SPECIFICATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axx = gca;
    
    ylim([0 70])
    ylabel('Altitude (km)')
    set(gca,'linewi',1.5); grid on;
    
    axx.GridColor = 'k';
    axx.TickDir = 'out';
    yminortick(0:5:100);
    
    switch ax
%         case 1
%             logscale('x')
%             xlabel('GWMF (mPa)')
%             set(gca,'xminorgrid','off')
%             xlim([0.01 1000])
%             
%             % draw xticks as real numbers
%             axx.XTick = [0.01 0.1 1 10 100 1000]; 
%             axx.XTickLabel = cellstr(num2str(axx.XTick(:)))';
            
        case 1
            
            xlim([0 6])
            axx.XTick = 0:10;
            xminortick(0:0.5:10)
            xlabel('Wave Amplitude (K)')
            
        case 2
            
            logscale('x')
            set(gca,'xgrid','off');
            xlabel('Wave Amplitude (K)')
            
%             xlim([0.06 25])
            xlim([0.1 10])
            
%             loglines = 10.^[-5:0.25:5];
            minorloglines = [0.002:0.001:0.009 0.02:0.01:0.09 0.2:0.1:0.9 2:1:9];
            majorloglines = [0.001 0.01 0.1 1 10];
            
            for i = 1:length(minorloglines)
                hold on; plot(minorloglines(i).*exp((zz-0)./(2*7)),zz,'color',[0 0 0 0.15],'linewi',2,'linest',':');
            end
            for i = 1:length(majorloglines)
                hold on; plot(majorloglines(i).*exp((zz-0)./(2*7)),zz,'color',[0 0 0 0.15],'linewi',2,'linest','-');
            end
            
            % draw xticks as real numbers
            axx.XTick = [0.1 0.2 0.5 1 2 5 10 20];
            axx.XTickLabel = cellstr(num2str(axx.XTick(:)))';
            
    end
    
    %%%% MODEL SPONGE LAYER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xl = get(gca,'xlim');
    hold on; plot(xl,[58.5 58.5],'color',rgbtrip(.2),'linest','--','linewi',1.5);
%     switch ax
%         case 2
%             text(1.1*xl(1),58.5 + 5.5,'Model Sponge Layer','fontsize',0.75*fs,'color',mcolor('blue'));
%     end
    
    % PLOT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     switch ax
%         case 1
%             hold on; PLOTS(3) = plot(smoothn(1000.*MA.MFv,3),zlev,'color',MAcol,'linewi',lw);
%             hold on; PLOTS(2) = plot(smoothn(1000.*M.MFv,3),zlev,'color',Mcol,'linewi',lw);
%             hold on; PLOTS(1) = plot(1000.*A.MFv(inds),zlev(inds),'color',Acol,'linewi',lw);
%         case {2,3}
%             hold on; PLOTS(3) = plot(smoothn(MA.Av,3),zlev,'color',MAcol,'linewi',lw);
%             hold on; PLOTS(2) = plot(smoothn(M.Av,3),zlev,'color',Mcol,'linewi',lw);
%             hold on; PLOTS(1) = plot(A.Av(inds),zlev(inds),'color',Acol,'linewi',lw);
%     end
%     switch ax
%         case 1
%             hold on; PLOTS(3) = plot(MA.(R).MFvec,zlev,'color',MAcol,'linewi',lw);
%             hold on; PLOTS(2) = plot( M.(R).MFvec,zlev,'color',Mcol,'linewi',lw);
%             hold on; PLOTS(1) = plot( A.(R).MFvec,zlev,'color',Acol,'linewi',lw);
%         case {2,3}
            hold on; PLOTS(3) = plot(MA.(R).Avec,zlev,'color',MAcol,'linewi',lw);
            hold on; PLOTS(1) = plot( M.(R).Avec,zlev,'color',Mcol,'linewi',lw);
            hold on; PLOTS(2) = plot( A.(R).Avec,zlev,'color',Acol,'linewi',lw);
%     end
    
    %%%% AIRS VERTICAL WINDOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch ax
        case 2
            axx = gca;
            axx2 = axes;
%             xlim([0.06 25])
            set(axx2,'color','none','position',axx.Position,'xlim',axx.XLim,'ylim',axx.YLim,'xcolor','none','ycolor','none','xtick',[],'ytick',[],'box','off')
            winvec(zlev > 70 | zlev <= 0) = NaN;
            xpoint = 10.25;
            hold on; plot(axx2,xpoint+(0.5*winvec),zlev,'color',mcolor('red'),'linewi',1.5);
            set(axx2,'clipping','off')
            text(xpoint+0.9,58,'AIRS Measurement Window','fontsize',0.8*fs,'rotation',-90,'color',mcolor('red'));
            axes(axx);
    end
    
    
    % LEGEND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L = legend(PLOTS,{'Model','AIRS','Model as AIRS'});
    L.Location = 'southeast';
    L.FontSize = 0.8*fs;
    
    % LABEL
    nph_text([0 1.02],['(' alphabet(ax) ')'],'fontsize',1.75*fs);
    
    
    setfont(fs)
    
end







return


%% EXPORT? ================================================================

savename = ['~/Desktop/AirsModel_ExpAmp'];

disp(['Exporting to ' savename '...'])

nph_saveas(gcf,savename,'png')

disp('Done.')







return




































































% TBP = {...
%     A.MFh,...
%     M.MFh,...
%     MA.MFh,...
%     A.MFh-MA.MFh,...
%     A.MFv,...
%     M.MFv,...
%     MA.MFv,...
%     A.MFv-MA.MFv};

TBP = {...
    A.Ginih,...
    M.Ginih,...
    MA.Ginih,...
    A.Ginih-MA.Ginih};

% Trim to vertical domains:
% TBP{1}(zlev < 18 | zlev >= 62,:) = NaN;
% TBP{2}(zlev < 7,:) = NaN;
% TBP{3}(zlev < 7,:) = NaN;

% % Trim to horizontal domains:
% xinds = abs(X') > 500;
% TBP{1}(xinds) = NaN;
% TBP{2}(xinds) = NaN;
% TBP{3}(xinds) = NaN;

% AIRS vertical window, as applied in 3DST:
% lower_height_index = 12; % corresponds to 12*1.5 = 18km
lower_height_index = 13; % corresponds to 13*1.5 = 19.5km
upper_height_index = 40; % corresponds to 40*1.5 = 60km
tapering_length = 6; % corresponds to 6*1.5 = 9km window at edge
OUT = sgwex_airs_apply_height_window(ones(60,80,50),lower_height_index,upper_height_index,tapering_length);
winvec = squeeze(OUT(1,1,:));
window = 1 - squeeze(OUT(30,:,:));
% window(xinds) = NaN;


% Specs:

cl = [0.2 0.8];

clims = {...
    cl,...
    cl,...
    cl,...
    cl};

xlims = {...
    [-500 500],...
    [-500 500],...
    [-500 500],...
    [-500 500]};

ylims = {...
    [-400 400],...
    [-400 400],...
    [-400 400],...
    [-400 400]};

tits = {...
    'AIRS',...
    'Model',...
    'Model as AIRS',...
    ['Difference' newline '(AIRS - Model as AIRS)']};

ct = 0:0.2:1;

cticks = {...
    ct,...
    ct,...
    ct,...
    ct};



%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % At the end, a vertical profile of Gini over the island:
% subplot(rows,cols,4)
%
% TBP = {A.Giniv,M.Giniv,MA.Giniv};
%
% % Trim AIRS to vertical domain
% TBP{1}(zlev < 20 | zlev > 60) = NaN;
%
% for i = 1:3
%     hold on; plot(TBP{i},zlev); grid on;
% end
%
% return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%



for ax = 1:3
    
    subplot(rows,cols,ax)
    
    tbp = TBP{ax}';
    
    switch ax
        case {1,2,3,4}
            [X,Y] = ndgrid(linspace(-600,600-15,80),linspace(-450,450,60));
        case {5,6,7,8}
            [X,Y] = ndgrid(linspace(-600,600-15,80),linspace(0,75,50));
            % Trim to region:
            % tbp(X < xlims{ax}(1) | X > xlims{ax}(2)+1) = NaN;
            % tbp(Y < ylims{ax}(1) | Y > ylims{ax}(2)+1) = NaN;
    end
    
    % Contour levels:
    %     switch ax
    %         case {2,6}
    %             ginilev = clims{ax}(1):5:clims{ax}(2);
    %         case {4,8}
    %             ginilev = (3*clims{ax}(1)):0.5:(3*clims{ax}(2));
    %         otherwise
    ginilev = clims{ax}(1):0.05:clims{ax}(2);
    %     end
    
    hold on; contourf(X,Y,tbp,ginilev,'edgecolor','none'); grid on;
    
    
    % COLOURS
    cmap = colormap(gca,flipud(cbrew('Blues',length(ginilev))));
    clim(clims{ax})
    
    
    xlim(xlims{ax})
    ylim(ylims{ax})
    
    axx = gca;
    axx.XTick = -600:200:600;
    axx.XMinorTick = 'on';
    axx.XAxis.MinorTickValues = -600:100:600;
    
    axx.YAxis.TickDirection = 'out';
    
    
    set(gca,'box','on','linewi',1.5,'layer','top')
    
    
    
    % % %     % AIRS VERTICAL WINDOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % %     switch ax
    % % %         case {5,7,8}
    % % %             % AIRS VERTICAL WINDOW:
    % % %             mask = zeros(size(tbp));
    % % %             hold on; p = pcolor(X,Y,mask);
    % % %             p.FaceColor = rgbtrip(.75);
    % % %             p.EdgeColor = 'none';
    % % %             p.AlphaData = window;
    % % %             p.FaceAlpha = 'flat';
    % % %         case {2,3,5,6}
    % % %             axx.GridColor = 'w';
    % % %     end
    
    
    % Add topography
    % % %     switch ax
    % % %
    % % %         case {1,2,3,4}
    % South Georgia
    LonBox = [-40 -30];
    LatBox = [-60 -30];
    C = nph_draw_coastline([LonBox(1) LatBox(1) ; LonBox(2) LatBox(2)],0,0,'noplot','color','k');
    for i = 1:length(C)
        if length(C(i).Lon) > 1
            [d,az] = distance(clat,clon,C(i).Lat,C(i).Lon);
            hold on; plot(deg2km(d).*sind(az),deg2km(d).*cosd(az),'color','w','linewi',2.5);
            hold on; plot(deg2km(d).*sind(az),deg2km(d).*cosd(az),'color','k','linewi',1.5);
        end
    end
    
    set(gca,'xtick',-400:200:400);
    set(gca,'ytick',-400:200:400);
    
    
    tit = nph_text([0.45 1.025],tits{ax},'fontsize',1.3*fs,'fontweight','bold','horizontalalignment','center');
    %             tit.VerticalAlignment = 'top';
    %             tit.Position = [1 1.075 1 1] .* get(gca,'Position');
    % % %             if ax == 4
    % % %                 tit.Position = [1 1.04 1 1] .* tit.Position;
    % % %             end
    
    
    % COLORBARS
    apos = get(gca,'Position');
    cbar = nph_colorbar;
    cbar.Location = 'southoutside';
    drawnow;
    set(gca,'Position',apos);
    cpos = cbar.Position .* [1 1 0.8 1.25];
    cbar.Label.String = 'Gini Coefficient';
    cbar.Ticks = cticks{ax};
    
    % % %         case {5,6,7,8}
    % % %
    % % %             axes(axx)
    % % %             xd = xdist([1:end 1]);
    % % %             tp = topo([1:end 1])./1000;
    % % %             hold on; fill(xd,tp,'k','edgecolor',[.9 .9 .9],'linewidth',1.5);
    % % %             hold on; plot(xlims{ax},[0 0],'color','k','linewi',1)
    % % %
    % % %             set(gca,'xtick',-400:200:400);
    % % %             set(gca,'ztick',0:10:70);
    % % %
    % % %     end
    
    switch ax
        
        case {1,2,3,4}
            ylabel('Altitude (km)')
            xlabel('Zonal Distance (km)')
            ytix = get(gca,'yticklabel');
            ytix{2} = [' ' ytix{2}];
            set(gca,'yticklabel',ytix);
            
        case {5,6,7,8}
            xlabel('Zonal Distance (km)')
            ylabel('Meridional Distance (km)')
    end
    
    
    %%%% AIRS VERTICAL WINDOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch ax
        case 3
            axx = gca;
            axx2 = axes;
            set(axx2,'color','none','position',axx.Position,'xlim',axx.XLim,'ylim',axx.YLim,'xcolor','none','ycolor','none','xtick',[],'ytick',[],'box','off')
            winvec(zlev > 70 | zlev <= 0) = NaN;
            xpoint = 6;
            hold on; plot(axx2,xpoint+(50*winvec),zlev,'color','k','linewi',1.5);
            set(axx2,'clipping','off')
            text(xpoint+1,56.5,'AIRS Measurement Window','fontsize',0.8*fs,'rotation',-90);
    end
    
    %     switch ax
    %         case {1,4}
    %             hold on; nph_text([0.075 0.925],['(' alphabet(ax) ')'],'fontsize',1.75*fs,'color','k','textborder','w');
    %         case {2,3,5,6}
    hold on; nph_text([0.075 0.925],['(' alphabet(ax) ')'],'fontsize',1.75*fs,'color','w','textborder','k');
    %     end
    
    
    
    setfont(fs)
    
    
    
    
end












return


%% EXPORT? ================================================================

savename = ['~/Desktop/AirsModel_Gini'];

disp(['Exporting to ' savename '...'])

nph_saveas(gcf,savename,'png')

disp('Done.')







return








return



















