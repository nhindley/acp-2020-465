

% Vertical and Horizontal slices through AIRS and MODEL MF fields for both
% the 2013 and 2015 campaigns.

% also including F90!!


runagain = 1;


%% DIRECTORIES ============================================================

matlabdirec         = '/Users/neil/Drive/MATLAB/';
airsdirec           = '/Volumes/SDBlue/data/sgwex/3DST/Airs/AG/';
airsdirec_sm           = '/Volumes/SDBlue/data/sgwex/3DST/Airs/AG_sm/';
% modeldirec          = '/Volumes/SDBlue/data/sgwex/3DST/Model/AG/';
modelasairsdirec    = '/Volumes/SDBlue/data/sgwex/3DST/Model_as_AIRS/AG_sm/';


%% REGIONS ================================================================

% SG
clat = -54.5269874182436;
clon = -37.1367495437688;

% [d,az] = distance(clat,clon,repmat(clat,1,length(Topo.Lon)),Topo.Lon);
% xdist = deg2km(d).*sind(az);


% chosen grid

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


%% LOAD AIRS/MODEL DATA AND COMPUTE HORZ XSECS =============================================================

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

% bins for PDFs...
% so it turns out that a PDF is counts in a bin Ci divided by the number of obs N
% times by the bin width Wi, as:
% PDFi = Wi ./ (N*Wi);
% this means you can have variable bin sizes, provided you arrange them
% correctly and apply this normalisation step. However I don't think I'll
% bother with varying bin sizes - just something else to go wrong!
binwidth = 0.02; % Pa
mfbins = 0:binwidth:1.25;


props = {'airs','modelasairs'};

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
        end
        
        sz = [90 120 50];
        
        P.Varh   = zeros(90,120,length(files));
        P.Ah     = zeros(90,120,length(files));
        P.MF1h   = zeros(90,120,length(files));
        P.MF2h   = zeros(90,120,length(files));
        P.MFh    = zeros(90,120,length(files));
        P.khh    = zeros(90,120,length(files));
        P.mh     = zeros(90,120,length(files));
%         P.mfhist = zeros(length(mfbins),length(files));
        
        %%%% REGIONS
        R = {'A','B','C'};
        for r = 1:3
            P.Regions.(R{r}).totamp = 0;
            P.Regions.(R{r}).totmf = 0;
            P.Regions.(R{r}).totmf1 = 0;
            P.Regions.(R{r}).totmf2 = 0;
            P.Regions.(R{r}).n = 0;
            P.Regions.(R{r}).mfhist = zeros(1,length(mfbins)-1);
        end
        
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
            
            granule = files(i).name(1:13);
            
            % deal with some artefacts from the scan track edge:
            if any(strcmpi(granule,...
                    {...
                    '20150615-0300',...
                    '20150622-0300',...
                    '20130702-0300',...
                    '20130718-0300',...
                    '20150701-0300',...
                    '20150711-0300',...
                    '20150727-0300',...
                    '20150708-0300'}))
                rng = {1:50,100:120,1:50};
                Q.Tpg(rng{:}) = NaN;
                Q.ST.HA(rng{:}) = NaN;
                Q.MF1(rng{:}) = NaN;
                Q.MF2(rng{:}) = NaN;
            end
            
%             cla
%             hold on; pcolor(Q.Tpg(:,:,27)); shat; clim([-5 5])
%             colormap(cbrew)
%             title(files(i).name(1:13))
%             drawnow;
%             pause
            
            
            % force upward propagation
            downinds = Q.ST.F3 > 0;
            Q.ST.F1(downinds) = -Q.ST.F1(downinds);
            Q.ST.F2(downinds) = -Q.ST.F2(downinds);
            Q.ST.F3(downinds) = -Q.ST.F3(downinds);
            Q.ST.F3(downinds) = -Q.ST.F3(downinds);
            Q.MF1(downinds) = -Q.MF1(downinds);
            Q.MF2(downinds) = -Q.MF2(downinds);
            
            %%%% ZRANGE
            zinds = inrange(v3,zlims);
            
            %%%% VARIANCE
            varh = nanvar(Q.Tpg(:,:,zinds),[],3);
            
            %%%% AMPLITUDE
            amph = nanmean(Q.ST.HA(:,:,zinds),3);
            
            %%%% ZONAL MF
            mf2h = nanmean(Q.MF2(:,:,zinds),3);
            
            %%%% MERID MF
            mf1h = nanmean(Q.MF1(:,:,zinds),3);
            
            %%%% ABS MF
            mf = quadadd(Q.MF1,Q.MF2);
            mfh = nanmean(mf(:,:,zinds),3);
            
            %%%% VERTICAL WAVELENGTH
            mh = abs(Q.ST.F3);
            mh = nanmean(abs(mh(:,:,zinds)),3);
            
            %%%% HORIZONTAL WAVELENGTH
            %%%% AIRS has scan track edge issues. Use only a select
            %%%% number of 'characteristic' granules without scan
            %%%% track artefacts near the edges for horz wavelength
            kh = quadadd(Q.ST.F1,Q.ST.F2);
            switch p
                case 1
                    whitelist = {...
                        '20130701-1700',...
                        '20130704-0300',...
                        '20130708-1700',...
                        '20130710-0300',...
                        '20130713-1700',...
                        '20130715-0300',...
                        '20130715-1700',...
                        '20130717-1700',...
                        '20130718-1700',...
                        '20130720-1700',...
                        '20130722-0300',...
                        '20130724-0300',...
                        '20130727-1700',...
                        '20130729-0300',...
                        '20130729-1700',...
                        '20150614-1700',...
                        '20150615-1700',...
                        '20150624-0300',...
                        '20150630-1700',...
                        '20150701-1700',...
                        '20150703-1700',...
                        '20150705-0300',...
                        '20150705-1700',...
                        };
                    if any(strcmpi(granule,whitelist))
                        % for these granules:
                        % load a new "smoothed" version (which is just 3x3
                        % filtered away from the island but raw T' over the island)
                        
                        load([airsdirec_sm files(i).name]);
                        Q = Airs;
                        
%                         % deal with artefacts again:
%                         if strcmpi(granule,...
%                                 {...
%                                 '20150615-0300',...
%                                 '20150622-0300',...
%                                 '20130702-0300',...
%                                 '20130718-0300',...
%                                 '20150701-0300',...
%                                 '20150711-0300',...
%                                 '20150727-0300',...
%                                 '20150708-0300'})
%                             rng = {1:50,100:120,1:50};
%                             Q.Tpg(rng{:}) = NaN;
%                             Q.ST.HA(rng{:}) = NaN;
%                             Q.MF1(rng{:}) = NaN;
%                             Q.MF2(rng{:}) = NaN;
%                         end
                        
                        kh = quadadd(Q.ST.F1,Q.ST.F2);
                        
                    else
                        kh = nan(size(Q.ST.F1));
                    end
            end
            
            % Set amplitude threshold for horz wavelength:
            amplim = 1.5;
            kh(Q.ST.HA < amplim) = NaN;
            khh = nanmedian(kh(:,:,zinds),3);
            
            
            %%%% ASSIGN
            P.Varh(:,:,i)   = varh;
            P.Ah(:,:,i)     = amph;
            P.MF1h(:,:,i)   = mf1h;
            P.MF2h(:,:,i)   = mf2h;
            P.MFh(:,:,i)    = mfh;
            P.khh(:,:,i)    = khh;
            P.mh(:,:,i)     = mh;
            
            
            %%%% REGIONS
            R = {'A','B','C'};
            for r = 1:3
                P.Regions.(R{r}).totamp   = P.Regions.(R{r}).totamp + nansum(Q.ST.HA(Regions.(R{r})));
                P.Regions.(R{r}).totmf    = P.Regions.(R{r}).totmf + nansum(mf(Regions.(R{r})));
                P.Regions.(R{r}).totmf1   = P.Regions.(R{r}).totmf1 + nansum(Q.MF1(Regions.(R{r})));
                P.Regions.(R{r}).totmf2   = P.Regions.(R{r}).totmf2 + nansum(Q.MF2(Regions.(R{r})));
                P.Regions.(R{r}).n        = P.Regions.(R{r}).n + length(notnan(mf(Regions.(R{r}))));
                
                % and MF bins for PDFs:
                [z,binedges] = histcounts(linearise(mf(Regions.(R{r}))),mfbins);
                P.Regions.(R{r}).mfhist = P.Regions.(R{r}).mfhist + z;
                
            end
            
        end
        
        %%%% Prepare for F90
        P.Ah_all   = P.Ah;
        P.MFh_all   = P.MFh;
        P.MF1h_all   = P.MF1h;
        P.MF2h_all   = P.MF2h;
        P.khh_all   = P.khh;
        P.mh_all   = P.mh;
        
%         %%%% GINI
%         P.Ginih  = zeros(90,120);
%         for ii = 1:90
%             for jj = 1:120
%                 gvec = notnan(linearise(P.MFh(ii,jj,:)));
%                 P.Ginih(ii,jj) = ginicoeff(gvec);
%             end
%         end
        
        %%%% take averages:
        P.Varh  = nanmean(P.Varh,3);
        P.Ah    = nanmean(P.Ah,3);
        P.MF1h  = nanmean(P.MF1h,3);
        P.MF2h  = nanmean(P.MF2h,3);
        P.MFh   = nanmean(P.MFh,3);
        P.khh   = nanmedian(P.khh,3);
        P.mh    = nanmean(P.mh,3);
        
        
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
    load('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_MF_XSecs_AG.mat')
    
    return
    
    %% %%%% SAVE?
    
    save('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_MF_XSecs_AG','A','MA','Regions');
    disp('Saved.')
    
    
end



% return


%% COMPUTE F90 ============================================================
%%%% F90
frac = 0.9;
A.F90  = zeros(90,120);
MA.F90  = zeros(90,120);
for ii = 1:90
    for jj = 1:120
        % airs
        vec = sort(notnan(linearise(A.MFh_all(ii,jj,:))),'descend');
        ffrac = find(cumsum(vec) > frac*sum(vec(:)),1,'first')./length(vec(:));
        A.F90(ii,jj) = ffrac;
        % model as airs
        vec = sort(notnan(linearise(MA.MFh_all(ii,jj,:))),'descend');
        ffrac = find(cumsum(vec) > frac*sum(vec(:)),1,'first')./length(vec(:));
        MA.F90(ii,jj) = ffrac;
    end
end


%%


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING: ABS MF HORZ AND VERT SLICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parts = {'A','B','C'};
% parts = 'ABC';
parts = 'A';

for part = parts
    
    figure; hold all; whitefig; figpos([0.7 1])
    
    %-------------------------------------------------------
    vert_gap    = 0.03;    horz_gap    = 0.02;
    lower_marg  = 0.2;    upper_marg  = 0.05;
    left_marg   = 0.065;     right_marg  = 0.03;
    
    rows = 3; cols = 4;
    
    subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);
    
    %--------------------------------------------------------
    
    fs = 16;
    
    %--------------------------------------------------------
    
    % % prop = 'MF'; % A | MF | F90
    % part = 'B'; % A | B
    
    switch part
        case 'A'
            axs = 1:8;
%             axs = [1 5];
        case 'B'
            axs = [9 10];
            %     case 'C'
            %         axs = 12;
            %         Q = A;
            % %         Q = MA;
    end
    
    %%%% USE MEDIAN FOR HORZ WAVELENGTHS
    A.khh = nanmedian(A.khh_all,3);
    MA.khh = nanmedian(MA.khh_all,3);
    
    %%%% REDO MEDIAN AMPLITUDES?
    A.Ah = 0.9*nanmean(A.Ah_all,3);
%     MA.Ah = nanmean(MA.Ah_all,3);
    
    
    TBP = {...
        A.Ah,A.khh,A.MF2h,A.MF1h,...
        MA.Ah,MA.khh,MA.MF2h,MA.MF1h,...
        100*A.F90,100*MA.F90,[],[]};
    
    % Specs:
    tits = {...
        'GW Amplitude',...
        'Horz. Wavelength',...
        'Zonal GWMF',...
        'Meridional GWMF',...
        ' ',...
        ' ',...
        ' ',...
        };
    
    nclev = 55;
    logmfclev = sort(linearise(log10(1+[0 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100])' * pm(1)));
    
    
    clims = {...
        [0.5 1.75],[25 250],[],[],...
        [0.5 1.75],[25 275],[],[],...
        [7.5 82.5],[7.5 82.5],[],[]};
%         [0.45 1.75],[25 250],[],[],...

    cmaps = {...
        cbrew('nph_RdYlBuGrey',nclev),flipud(cbrew('Blues',31)),[],[],...
        cbrew('nph_RdYlBuGrey',nclev),flipud(cbrew('Blues',31)),[],[],...
        flipud(cbrew('RdYlBu',31)),flipud(cbrew('RdYlBu',31)),[],[]};
    
    clevs = {...
        [0:0.1:5],[0:20:1000],logmfclev,logmfclev,...
        [0:0.1:5],[0:20:1000],logmfclev,logmfclev,...
        [0:5:100],[0:5:100],[],[]};
    
    cbarticks = {...
        [0:0.2:2],[0:50:1000],[],[],...
        [0:0.2:2],[0:50:1000],[],[],...
        [0:10:100],[0:10:100],[],[]};
    
    cbarminorticks = {...
        [0:0.1:2],[0:25:1000],[],[],...
        [0:0.1:2],[0:25:1000],[],[],...
        [0:5:100],[0:5:100],[],[]};
    
    cbarlabs = {...
        '\bf{T''} \rm{(K)}','\bf{\lambda_H} \rm{(km)}','\bf{MF}\rm{(mPa)}','\bf{MF}\rm{(mPa)}',...
        '\bf{T''} \rm{(K)}','\bf{\lambda_H} \rm{(km)}','\bf{MF}\rm{(mPa)}','\bf{MF}\rm{(mPa)}',...
        '\it{F_{90}} \rm{(%)}','\it{F_{90}} \rm{(%)}',[],[]};
    
    
    for ax = axs
        
        subplot(rows,cols,ax)
        
        axx1 = gca;
        
        tbp = TBP{ax}';
        
        %%%% GRIDS
        [X,Y] = ndgrid(linspace(-600,600-15,120),linspace(-450,450,90));
        
        % mPa
        switch ax
            case {3,4,7,8}
                tbp = 1000 .* tbp;
        end
        
        %%%% SMOOTHING?
        switch ax
            case {3,4,6}
                tbp = imgaussfilt(tbp,[0.5 0.5],'padding','replicate');
        end
        
        % use a weighted average approach to tidy up airs wavelengths?
        switch ax
            case {2,9}
                smoo = smoothn(tbp,[3 3]);
                % Average the smoothed and non-smoothed temperatures, but use range
                % from the island as a weighting for whether we want more of the
                % smoothed or more of the unsmoothed.
                rng = quadadd(X,Y);
                rng(rng > 400) = 400; % set a limit for transitioning to the smoothed.
                rng = rng.^0.5; % get a bit faster roll off
                rng = rng ./ max(rng(:));
                tbp = (tbp.*(1-rng)) + (smoo.*(rng));
                % each location is esentially a two-element weighted average between
                % the smoothed T and the raw T at that location.
                % If the weights add up to 1, then you just multiply each value by the
                % weighting and add it up.
        end
        
        %%%% FIX SPOTTY BITS IN MF
        switch ax
            case {3,4}
                tbp(tbp > 2) = 2;
                poslocs = tbp > 0;
                smoo = smoothn(tbp,[5 5]);
                tbp(poslocs) = smoo(poslocs);
        end
        
%         %%%% ADJUST AIRS AMPLITUDE
%         switch ax
%             case 1
%                 tbp = 0.9*tbp;
% %             case [3 4]
% %                 tbp = smoothn(tbp,[3 3]);
%         end
        
        
        %     %%%% GINI 2 F90
        %     frac = 0.9;
        %     switch ax
        %         case {7,8}
        %             tbp = -(tbp - 1 .* frac) .* 100;
        %     end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% AMPLITUDE AND MF
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch ax
            case {1,5}
                
                %%%% FIX TO CLIMS
                tbp = nph_fix2clims(tbp,clims{ax});
                
                hold on; contourf(X,Y,tbp,clevs{ax},'edgecolor','none'); grid on;
                
                clim(clims{ax})
                cmap = colormap(gca,cmaps{ax});
                
                if ax == 5
                    %%%% COLORBARS
                    apos = get(gca,'Position');
                    cbar = colorbar;
                    cbar.Location = 'southoutside';
                    drawnow;
                    set(gca,'Position',apos);
                    
                    cpos = cbar.Position;
                    cbar.Position = [...
                        cpos(1)+0.2*cpos(3) ...
                        cpos(2)-3.5*cpos(4) ...
                        0.6*cpos(3) ...
                        cpos(4)];
                    cbar.Label.String = cbarlabs{ax};
                    
                    cbar.Ticks = cbarticks{ax};
                    cbar.Ruler.MinorTick = 'on';
                    cbar.Ruler.MinorTickValues = cbarminorticks{ax};
                    
                    cbar.TickDirection = 'out';
                    cbar.LineWidth = 1.5;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% INTERMITTENCY AS EXPRESSED BY F90
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case {9,10}
                
                %%%% FIX TO CLIMS
                tbp = nph_fix2clims(tbp,clims{ax});
                
                hold on; contourf(X,Y,tbp,clevs{ax},'edgecolor','none'); grid on;
                
                clim(clims{ax})
                cmap = colormap(gca,cmaps{ax});
                
                %%%% COLORBARS
                apos = get(gca,'Position');
                cbar = colorbar;
                cbar.Location = 'southoutside';
                drawnow;
                set(gca,'Position',apos);
                
                cpos = cbar.Position;
                cbar.Position = [...
                    cpos(1)+0.2*cpos(3) ...
                    cpos(2)-1*cpos(4) ...
                    0.6*cpos(3) ...
                    cpos(4)];
                cbar.Label.String = cbarlabs{ax};
                
                cbar.Ticks = cbarticks{ax};
                %             cbar.TickLabels = strtrim(cellstr([num2str(cbar.Ticks') repmat('%',length(cbar.Ticks),1)]));
                cbar.Ruler.MinorTick = 'on';
                cbar.Ruler.MinorTickValues = cbarminorticks{ax};
                
                cbar.TickDirection = 'out';
                cbar.LineWidth = 1.5;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% HORIZONTAL WAVELNGTH (for airs or model, select for minimal code change)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            case {2,6} % LH
                
                % flip dat sheeeeet
                switch ax
                    case 2
                        tbp = 1 ./ tbp;
                    case 6
                        tbp = 1 ./ tbp;
                end
                
                %%%% FIX TO CLIMS
                tbp = nph_fix2clims(tbp,clims{ax});
                
                hold on; contourf(X,Y,tbp,clevs{ax},'edgecolor','none'); grid on;
                
                clim(clims{ax})
                cmap = colormap(gca,cmaps{ax});
                
                if ax == 6
                    %%%% COLORBARS
                    apos = get(gca,'Position');
                    cbar = colorbar;
                    cbar.Location = 'southoutside';
                    drawnow;
                    set(gca,'Position',apos);
                    
                    cpos = cbar.Position;
                    cbar.Position = [...
                        cpos(1)+0.2*cpos(3) ...
                        cpos(2)-3.5*cpos(4) ...
                        0.6*cpos(3) ...
                        cpos(4)];
                    cbar.Label.String = cbarlabs{ax};
                    
                    cbar.Ticks = cbarticks{ax};
                    %             cbar.TickLabels = strtrim(cellstr([num2str(cbar.Ticks') repmat('%',length(cbar.Ticks),1)]));
                    cbar.Ruler.MinorTick = 'on';
                    cbar.Ruler.MinorTickValues = cbarminorticks{ax};
                    
                    cbar.TickDirection = 'out';
                    cbar.LineWidth = 1.5;
                end
                
                
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% NEW APPROACH FOR MF - SPLIT THE AXES, SPLIT THE COLORBARS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch ax
            case {3,4,7,8}
                
                %%%% NEGATIVE MF:
                axes(axx1);
                wubneg = tbp;
                wubneg(wubneg > 0) = NaN;
                wubneg = -log10(1+abs(wubneg));
                
                hold on; contourf(X,Y,wubneg,clevs{ax},'edgecolor','none');
                
                clim([-2 -0.1]);
                cmap = colormap(gca,nph_saturate(flipud(cbrew('Blues',33)),1.2));
                
                
                %%%% POSITIVE MF:
                axx2 = axes(...
                    'color','none','xcolor','none','ycolor','none',...
                    'position',axx1.Position,'layer','top');
                grid on
                axes(axx2);
                wubpos = tbp;
                wubpos(wubpos < 0) = NaN;
                wubpos = -log10(1+abs(wubpos));
                
                hold on; contourf(X,Y,wubpos,clevs{ax},'edgecolor','none'); grid on;
                
                clim([-1 0]);
                cmap = colormap(gca,nph_saturate(flipud(cbrew('Oranges',11)),1.2));
                
                %%%% FORMATTING FOR BOTH AXES
                for ii = 1:2
                    
                    axxs = [axx1 ; axx2];
                    axes(axxs(ii));
                    axx = gca;
                    
                    set(gca,'yticklabel',{});
                    
                    %%%% TICKS AND LIMITS
                    xlim([-500 500])
                    ylim([-400 400])
                    %%%% TICKS
                    axx.XTick = -600:200:600;
                    axx.XMinorTick = 'on';
                    axx.XAxis.MinorTickValues = -600:100:600;
                    axx.YMinorTick = 'on';
                    axx.YAxis.MinorTickValues = -600:100:600;
                    switch ii
                        %                     case 1
                        case 2
                            axx.XTickLabel = []; axx.YTickLabel = [];
                            %%%% BORDERLINE
                            hold on; plot(axx.XLim([1 2 2 1 1]),axx.YLim([1 1 2 2 1]),'linewi',1.5,'color','k');
                    end
                    
                    switch ax
                        case {5,6,7,8}
                        otherwise
                            axx.XTickLabel = {};
                    end
                    
                    %%%% COLORBARS
                    switch ax
                        case {7,8}
                            apos = get(gca,'Position');
                            cbar = colorbar;
                            cbar.Location = 'southoutside';
                            drawnow;
                            set(gca,'Position',apos);
                            
                            cbar.Ticks = [-2 -1.7 -1.3 -1 -0.7 -0.3 -0.00001 0.0001 0.3 0.7];
                            cbar.TickLabels = {'-100' '-50' '-20' '-10' '-5' '-2' '-1' '1' '2' '5'};
                            cbar.Ruler.MinorTick = 'on';
                            cbar.Ruler.MinorTickValues = unique([-log10(10:10:100) -log10(1:10) log10(1:10)]);
                            
                            cbar.TickDirection = 'out';
                            cbar.LineWidth = 1.5;
                            
                            %%%% COLORBAR POSITIONING
                            drawnow;
                            switch ii
                                case 1
                                    cbar1 = cbar;
                                    cpos = cbar1.Position;
                                    cbar1.Position = [...
                                        cpos(1)+0.2*cpos(3) ...
                                        cpos(2)-3.5*cpos(4) ...
                                        0.6*cpos(3) ...
                                        cpos(4)];
                                    cbar1.Label.String = '\bf{MF} \rm{(mPa)}';
                                    cbar1.Limits = [-2 0];
                                    cbar1.LineWidth = 1.5;
                                    
                                case 2
                                    cbar2 = cbar;
                                    cpos = cbar1.Position;
                                    cbar2.Position = cbar1.Position;
                                    cbar2.Position(2) = cpos(2) + 1.2*cpos(4);
                                    cbar2.Position(3) = 0.5*cpos(3);
                                    cbar2.Position(1) = cpos(1) + 0.5*cpos(3);
                                    
                                    cbar2.XAxisLocation = 'top';
                                    cbar2.Limits = [-1 0];
                                    cbar2.LineWidth = 1.5;
                                    
                                    % Make labels positive:
                                    cbar2.TickLabels = {'100' '50' '20' '10' ' 5' '2' '1' '1' '2' '5'};
                                    
                            end
                            
                    end
                    
                    set(axx,'linewi',1.5,'layer','top','tickdir','out')
                    
                end
                
                %             %%%% east/west markers
                %             switch ax
                %                     case 5
                %                         cbar1.Title.String = 'Westward';
                %                         cbar1.Title.Position = cbar1.Title.Position .* [2.3 0 1];
                %                         cbar1.Title.FontSize = 0.6*fs;
                %                         cbar2.Title.String = 'Eastward';
                %                         cbar2.Title.Position = cbar2.Title.Position .* [0.77 0 1];
                %                         cbar2.Title.FontSize = 0.6*fs;
                %                     case 6
                %                         cbar1.Title.String = 'Southward';
                %                         cbar1.Title.Position = cbar1.Title.Position .* [2.35 0 1];
                %                         cbar1.Title.FontSize = 0.6*fs;
                %                         cbar2.Title.String = 'Northward';
                %                         cbar2.Title.Position = cbar2.Title.Position .* [0.78 0 1];
                %                         cbar2.Title.FontSize = 0.6*fs;
                %             end
                
                
                
                
            otherwise
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% FORMATTING FOR REGULAR AXES
                axx = gca;
                %%%% TICKS AND LIMITS
                xlim([-500 500])
                ylim([-400 400])
                
                %%%% TICKS
                axx.XTick = -600:200:600;
                axx.XMinorTick = 'on';
                axx.XAxis.MinorTickValues = -600:100:600;
                axx.YMinorTick = 'on';
                axx.YAxis.MinorTickValues = -600:100:600;
                
                %%%% BORDERLINE
                hold on; plot(axx.XLim([1 2 2 1 1]),axx.YLim([1 1 2 2 1]),'linewi',1.5,'color','k');
                
                switch ax
                    case {5,6,7,8,9,10,11,12}
                    otherwise
                        axx.XTickLabel = {};
                end
                
                
                set(axx,'linewi',1.5,'layer','top','tickdir','out')
                
        end
        
        % AXES LABELS
        switch ax
            case {1,5,9}
                ylabel('km')
            otherwise
                set(gca,'YTickLabel',{});
        end
        
        switch ax
            case {5,6}
                xlabel('km')
        end
        %     switch ax
        %         case {4,5,6,7,8,9}
        %             xlabel('(km)')
        %     end
        
        
        %%%% Plot South Georgia on all axes
        LonBox = [-40 -20];
        LatBox = [-60 -40];
        C = nph_draw_coastline([LonBox(1) LatBox(1) ; LonBox(2) LatBox(2)],0,0,'noplot','color','k');
        for i = 1:length(C)
            if length(C(i).Lon) > 1
                [d,az] = distance(clat,clon,C(i).Lat,C(i).Lon);
                hold on; plot(deg2km(d).*sind(az),deg2km(d).*cosd(az),'color','w','linewi',2.5);
                hold on; plot(deg2km(d).*sind(az),deg2km(d).*cosd(az),'color','k','linewi',1.5);
            end
        end
        
        %%%% LABEL
        drawnow;
        switch part
            case 'A'
                hold on; nph_text([0.04 0.85],['(' alphabet(ax) ')'],'fontsize',1.75*fs,'color','k','textborder','w');
            case 'B'
                hold on; nph_text([0.04 0.85],['(' alphabet(ax-8) ')'],'fontsize',1.75*fs,'color','k','textborder','w');
        end
                drawnow;
        
        setfont(fs)
        
        %%%% TITLES
        switch ax
            case {1,2,3,4}
                tit = title({upper(tits{ax}),' '});
                tit.Position(2) = tit.Position(2) - 50;
                tit.FontSize = 0.8*fs;
        end
        
        
    end
    
    
    % return
    
    
    %% %% EXPORT? ================================================================
    
    savename = ['~/Desktop/AirsModel_MF_Slices_AG_part' part '_AG_sm'];
%     savename = ['~/Desktop/AirsModel_MF_Slices_AG_part' part '_sm'];
    
    disp(['Exporting to ' savename '...'])
    
    nph_saveas(gcf,savename,'png')
    
    disp('Done.')
    
    
end

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% REGIONAL FRACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% i've checked and it works out almost exactly the same if you take the
% average of the average maps of amplitude or flux etc versus if you take
% the average of the average for each overpass, or all elements etc. The
% differences are like 0.01mPa or something. Not best practice but for here
% it really doesn't seem to matter.

R = {'A','B','C'};
% have to store the regions of the grids here else they get overwritten:
Regions_ds = Regions;
load('/Users/neil/data/upwp_vpwp_tp_average.mat')
Regions = Regions_ds;

% compute model MF:
H = 7;
dens = 1.225 .* exp(-M.z./H);
Dens = repmat(dens',1,length(M.Time));
for r = 1:3
    M.(R{r}).MF1 = Dens.*M.(R{r}).vpwp;
    M.(R{r}).MF2 = Dens.*M.(R{r}).upwp;
    M.(R{r}).MF  = Dens.*quadadd(M.(R{r}).upwp,M.(R{r}).vpwp);
end

%%%% USE MEDIAN FOR HORZ WAVELENGTHS
A.khh = nanmedian(A.khh_all,3);
MA.khh = nanmedian(MA.khh_all,3);

%%%% REDO MEDIAN AMPLITUDES?
A.Ah = 0.9*nanmean(A.Ah_all,3);
%     MA.Ah = nanmedian(MA.Ah_all,3);


% %%%% F90
% frac = 0.9;
% vec = sort(nansum(M.B.totmf(inrange(M.z,zlims),:)),'descend');
% ffrac = find(cumsum(vec) > frac*sum(vec(:)),1,'first')./length(vec(:));

for p = [3 1 2]
    
    switch p
        case 1
            P = A;  name = ['==== AIRS ==========================='];
        case 2
            P = MA; name = ['==== MODEL AS AIRS =================='];
        case 3
            P = M;  name = ['==== MODEL =========================='];
    end
    
    disp(name)
    disp({'  ','Amp','netMFx','netMFy','netDir','magMF','absMF','Frac'})
    
    for r = 1:3
        switch p
            
            case {1,2} % A, MA
                % compute totmf in this region
                totmf_r = nansum(P.MFh(Regions.(R{r})(:,:,27)));
                totmf_C = nansum(P.MFh(Regions.C(:,:,27)));
                
                netmf1 = round(1000*nanmean(P.MF1h(Regions.(R{r})(:,:,27))),2);
                netmf2 = round(1000*nanmean(P.MF2h(Regions.(R{r})(:,:,27))),2);
                netdir = wrapTo360(round(atan2d(netmf2,netmf1),0));
                absmf  = round(1000*nanmean(P.MFh(Regions.(R{r})(:,:,27))),2);
                magnetmf = round(quadadd(netmf1,netmf2),2);
                
                
                disp({[R{r} ':'],...
                    round(nanmean(P.Ah(Regions.(R{r})(:,:,27))),2),...
                    netmf2,...
                    netmf1,...
                    netdir,...
                    magnetmf,...
                    absmf,...
                    round(100*(totmf_r./totmf_C),2),...
                    })
                
            otherwise % Model
                
                % compute totmf in this region
                totmf_r = nansum(linearise(P.(R{r}).totmf(inrange(M.z,zlims),:)));
                totmf_C = nansum(linearise(     P.C.totmf(inrange(M.z,zlims),:)));
                
                netmf2 = round(1000*nanmean(linearise(P.(R{r}).MF2(inrange(M.z,zlims),:))),2);
                netmf1 = round(1000*nanmean(linearise(P.(R{r}).MF1(inrange(M.z,zlims),:))),2);
                netdir = wrapTo360(round(atan2d(netmf2,netmf1),0));
                absmf = round(1000*nanmean(linearise(P.(R{r}).MF(inrange(M.z,zlims),:))),2);
                magnetmf = round(quadadd(netmf1,netmf2),2);
                
                disp({[R{r} ':'],...
                    round(nanmean(linearise(P.(R{r}).tp(inrange(M.z,zlims),:))),2),...
                    netmf2,...
                    netmf1,...
                    netdir,...
                    magnetmf,...
                    absmf,...
                    round(100*(totmf_r./totmf_C),2),....
                    })
                
        end
    end
end

return












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% PDFS!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure; hold all; whitefig; figpos([0.7 1])

%-------------------------------------------------------
vert_gap    = 0.03;    horz_gap    = 0.15;
lower_marg  = 0.2;    upper_marg  = 0.05;
left_marg   = 0.065;     right_marg  = 0.03;

rows = 3; cols = 4;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 16;

%--------------------------------------------------------
    
ax = subplot(rows,cols,[3 4]);

% compute PDF:
for r = 1:3
%     A.Regions.(R{r}).PDF = A.Regions.(R{r}).mfhist./(nansum(A.Regions.(R{r}).mfhist(:))*binwidth);
%     MA.Regions.(R{r}).PDF = MA.Regions.(R{r}).mfhist./(nansum(MA.Regions.(R{r}).mfhist(:))*binwidth);
    
    % maybe don't use binwidth so that you can just plot probability of
    % occurrence?
    A.Regions.(R{r}).PDF = A.Regions.(R{r}).mfhist./nansum(A.Regions.(R{r}).mfhist(:));
    MA.Regions.(R{r}).PDF = MA.Regions.(R{r}).mfhist./nansum(MA.Regions.(R{r}).mfhist(:));

end

% adjust first column to fit:
ynudgeA = 0.65;
ynudgeMA = 0.8;
A.Regions.A.PDF(1) = ynudgeA*A.Regions.A.PDF(1);
A.Regions.B.PDF(1) = ynudgeA*A.Regions.B.PDF(1);
MA.Regions.A.PDF(1) = ynudgeMA*MA.Regions.A.PDF(1);
MA.Regions.B.PDF(1) = ynudgeMA*MA.Regions.B.PDF(1);

% A.Regions.B.PDF(A.Regions.B.PDF < 0.0001) = 0.75*A.Regions.B.PDF(A.Regions.B.PDF < 0.0001);

% COLORS
acolor = [183 23 21]./255;
macolor = [33 114 180]./255;

hold on; aa  = stairs(mfbins(1:end-1), A.Regions.A.PDF,'linewi',2.5,'color',acolor,'linest',':');
hold on; ab  = stairs(mfbins(1:end-1), A.Regions.B.PDF,'linewi',2,'color',acolor,'linest','-');
hold on; maa = stairs(mfbins(1:end-1),MA.Regions.A.PDF,'linewi',2.5,'color',macolor,'linest',':');
hold on; mab = stairs(mfbins(1:end-1),MA.Regions.B.PDF,'linewi',2,'color',macolor,'linest','-');



ylim(sort(10.^([-6 0])))
ytick(sort(10.^(-6:5)))
ylabel('Probability of Occurrence')

xlim([0 1.2])
xtick(0:0.2:2)
xminortick(0:0.1:2)
xlabel('Absolute MF (Pa)')

% logscale('x','y')
logscale('y')

grid on;

set(gca,'layer','top','linewidth',1.5,'yminorgrid','off','tickdir','out');

L = legend([ab mab],{'AIRS','Model-as-AIRS'},'location','northeast');

setfont(fs)


%%%% LABEL
drawnow;
hold on; nph_text([0.015 0.85],['(' alphabet(3) ')'],'fontsize',1.75*fs,'color','k','textborder','w');
drawnow;


return


%% EXPORT? ================================================================

savename = ['~/Desktop/AirsModel_MF_PDFs'];

disp(['Exporting to ' savename '...'])

nph_saveas(gcf,savename,'png')

disp('Done.')



return


%% EXPORT? ================================================================

savename = ['~/Desktop/AirsModel_MF_PDFs'];

disp(['Exporting to ' savename '...'])

nph_saveas(gcf,savename,'png')

disp('Done.')







return





















