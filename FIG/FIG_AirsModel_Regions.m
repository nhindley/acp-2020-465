
% What are the properties of GWs upwind and downwind of the island?
% to do this, two equal area boxes are taken which have cunning shapes.

% The output is a figure showing the regions and the mean/total fluxes,
% amplitudes, gini coefficents in both boxes for the A, M and MA.


%% EQUAL AREA BOXES =======================================================

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

% % model domain:
% h = 825;
% r = h/2;
% L = 1100;
% model domain:
h = 400;
r = h/2;
L = 600;

d1 = (pi*h/8) - (h/4) + ((L-r)/2);
d2 = L - r - d1;

%%% CREATE GRIDS ===========================================================

[D1,D2,D3] = ndgrid(...
    linspace(-450,450-15,60),...
    linspace(-600,600-15,80),...
    1.5:1.5:75);

ring = quadadd(D1,D2-(d1+r-(L/2))) <= r;

zlims = [25 55];

% keep awayyy from the edges
edgebox = inrange(D1,pm(h/2)) & inrange(D2,pm(L/2));

% Region A - Upwind of the island
Regions.A = ~ring & edgebox & inrange(D2,[-(L/2) (L/2)-d2]) & inrange(D3,zlims);
% Regions.A = (~ring | D2 <= (d1-(L/2))) & ~(D2 >= ((L/2)-d2)) & abs(D1) <= h/2;
% Regions.A = Regions.A & inrange(D3,zlims);

% Region B - Over/Downwind of the island
Regions.B = (ring | inrange(D2,[(L/2)-d2 (L/2)])) & edgebox & inrange(D3,zlims);
% Regions.B = ring | inrange(D2,[(d1+r-(L/2)) (L/2)]) & abs(D1) <= h/2;
% Regions.B = Regions.B & inrange(D3,zlims);

Regions.C = Regions.A | Regions.B;


%% DIRECTORIES ============================================================

runagain = 1;

load('/Users/neil/data/20150705-1700_AirsSG_3DST.mat')
sfvec = Airs.AltAmpScaling;
sf = permute(repmat(sfvec,1,60,80),[2 3 1]);

% airsdirec           = '/Volumes/SDBlue/data/sgwex/3DST/Airs/downsized_new/';
% modeldirec          = '/Volumes/SDBlue/data/sgwex/3DST/Model/downsized_new/';
% modelasairsdirec    = '/Volumes/SDBlue/data/sgwex/3DST/Model_as_AIRS/downsized_new/';
airsdirec           = '/Volumes/SDBlue/data/sgwex/3DST/Airs/4dp/downsized/';
modeldirec          = '/Volumes/SDBlue/data/sgwex/3DST/Model/4dp/downsized/';
modelasairsdirec    = '/Volumes/SDBlue/data/sgwex/3DST/Model_as_AIRS/4dp/downsized/';



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

xbox = plusminus(650);
ybox = plusminus(500);

xvec = xbox(1):1.5:xbox(2);
yvec = ybox(1):2.5:ybox(2);

[Y,X] = ndgrid(yvec,xvec);

% Convert distance box to Lat Lon:
[LAT,LON] = reckon(clat,clon,km2deg(quadadd(X,Y)),atan2d(X,Y));

% pull down to region:
tpr = Ftp(LAT,LON) ./ 1000; tpr(tpr < 0) = 0;
Mapr = uint8(cat(3,FMap1(LAT,LON),FMap2(LAT,LON),FMap3(LAT,LON))) * 1.5; % brighten map

% % Grayscale?
% Mapr = rgb2gray(Mapr);

% Coastline?
boundingbox = [min(LON(:)) min(LAT(:)) ; max(LON(:)) max(LAT(:))];
C = nph_draw_coastline(boundingbox,0,'noplot');

% smooth topography?
tpr = smoothn(tpr,[3 3]);


return



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
    
    props = {...
        'a','mf','totmf','totmf1pos','totmf1neg','totmf2pos','totmf2neg','Time'};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% LOAD AIRSO, MODEL AND MODELO AS AIRSO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    plotflag = 0;
    
    if plotflag
        figure; hold all; figpos([1 0.6]); whitefig;
    end
    
    % for each dataset...
    for q = 1:3
        
        % ASSIGN OUTPUT STRUCTURE
        OUT = struct;
        
        % less code to assign fields to structs in a loop:
        for p = 1:length(props)
            
            pp = props{p};
            
            switch pp
                case {'a','mf','totmf','totmf1pos','totmf1neg'}
                    OUT.(pp).A = zeros(1,length(files{q}));
                    OUT.(pp).B = zeros(1,length(files{q}));
                    OUT.(pp).An = zeros(1,length(files{q}));
                    OUT.(pp).Bn = zeros(1,length(files{q}));
%                 case {'totmf','totmf1pos','totmf1neg'}
%                     OUT.(pp).A = 0;
%                     OUT.(pp).B = 0;
%                     OUT.(pp).An = 0;
%                     OUT.(pp).Bn = 0;
            end
            
        end
        
        files{q}
        
        % NOW LOAD FILES
        for i = 1:length(files{q})
            
            disp(files{q}(i).name)
            load([direcs{q} files{q}(i).name])
            
            switch q
                case 1
                    Q = Airs;
                case {2,3}
                    Q = Model;
            end
            
            
            if plotflag
                
%             subplot(1,2,1)
%             cla
%             hold on; pcolor(Q.ST.HR(:,:,27)); shat;
%             colormap(cbrew);
%             clim([-5 5])
%             axis square
%             title(Q.Time)
%             
%             subplot(1,2,2)
%             cla
%             hold on; pcolor(sq(Q.ST.HR(30,:,:))'); shat;
%             colormap(cbrew);
%             clim([-5 5])
%             axis square
%             title(Q.Time)

                        

            
            subplot(1,2,2)
            cla
%             try
%             load([direcs{3} files{q}(i).name(1:13) '_Model_as_AIRS_3DST_ds_4dp'])
%             catch err
%                 err
%                 continue
%             end
            Q.Tp = Q.Tpg_sc ./ permute(repmat(Q.AltAmpScaling,1,60,80),[2 3 1]);
            Q.bg = Q.Tg - Q.Tp;
            hold on; pcolor(Q.bg(:,:,27)); shat;
%             hold on; pcolor(Model.Tg(:,:,27)); shat;
%             colormap(cbrew);
            cl = get(gca,'clim');
            colorbar;
            axis square
            title(Q.Time)
            
            

            subplot(1,2,1)
            cla
            hold on; pcolor(Q.Tg(:,:,27)); shat;
%             colormap(cbrew);
            clim(cl)
            colorbar;
            axis square
            title(Q.Time)
            
            setfont(16)
            
            drawnow;
            
            

            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % try UPWARD assumption instead?
            posinds = Q.ST.F3 > 0;
            Q.ST.F1(posinds) = -Q.ST.F1(posinds);
            Q.ST.F2(posinds) = -Q.ST.F2(posinds);
            Q.ST.F3          = -abs(Q.ST.F3);
            Q.MF1(posinds)   = -Q.MF1(posinds);
            Q.MF2(posinds)   = -Q.MF2(posinds);
            
            W = struct;
            %%%% Amplitude:            
            W.a.A   = nansum(linearise(Q.ST.HA(Regions.A) ./ sf(Regions.A)));
            W.a.B   = nansum(linearise(Q.ST.HA(Regions.B) ./ sf(Regions.B)));
            W.a.C   = nansum(linearise(Q.ST.HA(Regions.C) ./ sf(Regions.C)));
            W.a.An  = nansum(~isnan(Q.ST.HA(Regions.A)));
            W.a.Bn  = nansum(~isnan(Q.ST.HA(Regions.B)));
            W.a.Cn  = nansum(~isnan(Q.ST.HA(Regions.C)));
            
            %%%% Fluxes:
            amplim = 0.5;
            inds = ~isnan(Q.Tpg_sc) & Q.ST.HA >= amplim;
            
            % mean abs flux
            mf = quadadd(Q.MF1,Q.MF2);
            
%             W.mf.A = nanmean(linearise(mf(Regions.A & inds)));
%             W.mf.B = nanmean(linearise(mf(Regions.B & inds)));
%             W.mf.C = nanmean(linearise(mf(Regions.C & inds)));
%             W.mf.An  = nansum(~isnan(mf(Regions.A & inds)));
%             W.mf.Bn  = nansum(~isnan(mf(Regions.B & inds)));
%             W.mf.Cn  = nansum(~isnan(mf(Regions.C & inds)));
%             
%             % total fluxes
%             W.totmf.A = nansum(linearise(mf(Regions.A & inds)));
%             W.totmf.B = nansum(linearise(mf(Regions.B & inds)));
%             W.totmf.C = nansum(linearise(mf(Regions.C & inds)));
%             W.totmf.An = nansum(double(linearise(Regions.A & inds)));
%             W.totmf.Bn = nansum(double(linearise(Regions.B & inds)));
%             W.totmf.Cn = nansum(double(linearise(Regions.C & inds)));
%             
%             % mean merid flux
%             W.totmf1pos.A = nansum(Q.MF1(Q.MF1 > 0 & Regions.A & inds));
%             W.totmf1neg.A = nansum(Q.MF1(Q.MF1 < 0 & Regions.A & inds));
%             W.totmf1pos.B = nansum(Q.MF1(Q.MF1 > 0 & Regions.B & inds));
%             W.totmf1neg.B = nansum(Q.MF1(Q.MF1 < 0 & Regions.B & inds));
%             W.totmf1pos.C = nansum(Q.MF1(Q.MF1 > 0 & Regions.C & inds));
%             W.totmf1neg.C = nansum(Q.MF1(Q.MF1 < 0 & Regions.C & inds));
%             
%             W.totmf1pos.An = sum(linearise(Q.MF1 > 0 & Regions.A & inds));
%             W.totmf1neg.An = sum(linearise(Q.MF1 < 0 & Regions.A & inds));
%             W.totmf1pos.Bn = sum(linearise(Q.MF1 > 0 & Regions.B & inds));
%             W.totmf1neg.Bn = sum(linearise(Q.MF1 < 0 & Regions.B & inds));
%             W.totmf1pos.Cn = sum(linearise(Q.MF1 > 0 & Regions.C & inds));
%             W.totmf1neg.Cn = sum(linearise(Q.MF1 < 0 & Regions.C & inds));
%             
            % mean northward/southward flux
            ABC = {'A','B','C'};
            for d = 1:3
                abc = ABC{d};
                
                % mean abs flux (for GINI)
                W.mf.(abc) = nanmean(linearise(mf(Regions.(abc) & inds)));
                W.mf.([abc 'n'])  = nansum(~isnan(mf(Regions.(abc) & inds)));
                
                % total fluxes:
                W.totmf.(abc) = nansum(linearise(mf(Regions.(abc) & inds)));
                W.totmf.([abc 'n']) = nansum(double(linearise(Regions.(abc) & inds)));
                
                % total meridional fluxes:
                W.totmf1pos.(abc) = nansum(Q.MF1(Q.MF1 > 0 & Regions.(abc) & inds));
                W.totmf1neg.(abc) = nansum(Q.MF1(Q.MF1 < 0 & Regions.(abc) & inds));
                W.totmf1pos.([abc 'n']) = sum(linearise(Q.MF1 > 0 & Regions.(abc) & inds));
                W.totmf1neg.([abc 'n']) = sum(linearise(Q.MF1 < 0 & Regions.(abc) & inds));
                
                % total zonal fluxes:
                W.totmf2pos.(abc) = nansum(Q.MF2(Q.MF2 > 0 & Regions.(abc) & inds));
                W.totmf2neg.(abc) = nansum(Q.MF2(Q.MF2 < 0 & Regions.(abc) & inds));
                W.totmf2pos.([abc 'n']) = sum(linearise(Q.MF2 > 0 & Regions.(abc) & inds));
                W.totmf2neg.([abc 'n']) = sum(linearise(Q.MF2 < 0 & Regions.(abc) & inds));

            end
            
            
            %%%% ASSIGN
            
            W.Time = datenum(Q.Time);
            
            for p = 1:length(props)
                
                pp = props{p};
                
                switch pp
                    
%                     case {'a','mf'}
%                         
%                         OUT.(pp).A(i) = W.(pp).A;
%                         OUT.(pp).B(i) = W.(pp).B;
                        
                    case {'a','mf','totmf','totmf1pos','totmf1neg','totmf2pos','totmf2neg'}
                        
                        OUT.(pp).A(i) = W.(pp).A;
                        OUT.(pp).B(i) = W.(pp).B;
                        OUT.(pp).C(i) = W.(pp).C;
                        
                        OUT.(pp).An(i) = W.(pp).An;
                        OUT.(pp).Bn(i) = W.(pp).Bn;
                        OUT.(pp).Cn(i) = W.(pp).Cn;
                        
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
    
    
    
    % or load saved version
else
    
%     load('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_Regions')
    load('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_Regions_4dp')
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY FLUX FRACTIONS ETC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A.fractotmf.A = nansum(A.totmf.A) ./ nansum(A.totmf.A + A.totmf.B);
A.fractotmf.B = nansum(A.totmf.B) ./ nansum(A.totmf.A + A.totmf.B);
A.fractotmf.C = nansum(A.totmf.C) ./ nansum(A.totmf.C);
M.fractotmf.A = nansum(M.totmf.A) ./ nansum(M.totmf.A + M.totmf.B);
M.fractotmf.B = nansum(M.totmf.B) ./ nansum(M.totmf.A + M.totmf.B);
M.fractotmf.C = nansum(M.totmf.C) ./ nansum(M.totmf.C);
MA.fractotmf.A = nansum(MA.totmf.A) ./ nansum(MA.totmf.A + MA.totmf.B);
MA.fractotmf.B = nansum(MA.totmf.B) ./ nansum(MA.totmf.A + MA.totmf.B);
MA.fractotmf.C = nansum(MA.totmf.C) ./ nansum(MA.totmf.C);
A.gini.A = ginicoeff(A.mf.A);
A.gini.B = ginicoeff(A.mf.B);
A.gini.C = ginicoeff(A.mf.C);
M.gini.A = ginicoeff(M.mf.A);
M.gini.B = ginicoeff(M.mf.B);
M.gini.C = ginicoeff(M.mf.C);
MA.gini.A = ginicoeff(MA.mf.A);
MA.gini.B = ginicoeff(MA.mf.B);
MA.gini.C = ginicoeff(MA.mf.C);

names = {'AIRS','Model','Model as AIRS'};
data = {A,M,MA};

% for w = 1:3
%     W = data{w};
%     
%     disp(['%%%%%%%%%%%' names{w} '%%%%%%%%%%%%%%%%%%%%%%%%%'])
%     disp({'Mean Amp',round(nansum(W.a.A)./nansum(W.a.An),2),round(nansum(W.a.B)./nansum(W.a.Bn),2),round(nansum(W.a.C)./nansum(W.a.Cn),2)})
%     disp('     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%     disp({'Eastward Flux', round(1000*W.totmf2pos.A / (W.totmf2pos.An + W.totmf2neg.An),2),round(1000*W.totmf2pos.B / (W.totmf2pos.Bn + W.totmf2neg.Bn),2),round(1000*(W.totmf2pos.A + W.totmf2neg.A + W.totmf2pos.B + W.totmf2neg.B) / (W.totmf2pos.An + W.totmf2neg.An + W.totmf2pos.Bn + W.totmf2neg.Bn),2)})
%     disp({'Westward Flux', round(1000*W.totmf2neg.A / (W.totmf2pos.An + W.totmf2neg.An),2),round(1000*W.totmf2neg.B / (W.totmf2pos.Bn + W.totmf2neg.Bn),2),round(1000*(W.totmf2pos.A + W.totmf2neg.A + W.totmf2pos.B + W.totmf2neg.B) / (W.totmf2pos.An + W.totmf2neg.An + W.totmf2pos.Bn + W.totmf2neg.Bn),2)})
%     disp({'Northward Flux',round(1000*W.totmf1pos.A / (W.totmf1pos.An + W.totmf1neg.An),2),round(1000*W.totmf1pos.B / (W.totmf1pos.Bn + W.totmf1neg.Bn),2),round(1000*(W.totmf1pos.A + W.totmf1neg.A + W.totmf1pos.B + W.totmf1neg.B) / (W.totmf1pos.An + W.totmf1neg.An + W.totmf1pos.Bn + W.totmf1neg.Bn),2)})
%     disp({'Southward Flux',round(1000*W.totmf1neg.A / (W.totmf1pos.An + W.totmf1neg.An),2),round(1000*W.totmf1neg.B / (W.totmf1pos.Bn + W.totmf1neg.Bn),2),round(1000*(W.totmf1pos.A + W.totmf1neg.A + W.totmf1pos.B + W.totmf1neg.B) / (W.totmf1pos.An + W.totmf1neg.An + W.totmf1pos.Bn + W.totmf1neg.Bn),2)})
%     disp('     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%     disp({'Mean Flux',round(1000*nansum(W.totmf.A)./nansum(W.totmf.An),2),round(1000*nansum(W.totmf.B)./nansum(W.totmf.Bn),2),round(1000*nansum(W.totmf.C)./nansum(W.totmf.Cn),2)})
%     disp({'Total Flux %',round(W.fractotmf.A,2),round(W.fractotmf.B,2),[]})
%     disp({'Gini Coeff',round(W.gini.A,2),round(W.gini.B,2),round(W.gini.C,2)})
%     
% end

%% TRY AGAIN?

R = {'A','B','C'};
Rn = {'An','Bn','Cn'};

% for each dataset
for w = 1:3
    
    W = data{w};
    
    meanamp = struct;
    eastflux = struct;
    westflux = struct;
    zonalnet = struct;
    northflux = struct;
    southflux = struct;
    meridnet = struct;
    fractot = struct;
    
    % for each region
    for r = 1:3
        
        reg = R{r};
        regn = Rn{r};
    
    
    
    %     disp('     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    meanamp.(reg) = nansum(W.a.(reg))./nansum(W.a.(regn));
    
    num_in_region = nansum(W.totmf2pos.(regn) + W.totmf2neg.(regn));
    eastflux.(reg) = 1000 * (nansum(W.totmf2pos.(reg)) ./ num_in_region);
    westflux.(reg) = 1000 * (nansum(W.totmf2neg.(reg)) ./ num_in_region);
    zonalnet.(reg) = 1000 * (nansum(W.totmf2pos.(reg) + W.totmf2neg.(reg)) ./ num_in_region);
    
    num_in_region = nansum(W.totmf1pos.(regn) + W.totmf1neg.(regn));
    northflux.(reg) = 1000 * (nansum(W.totmf1pos.(reg)) ./ num_in_region);
    southflux.(reg) = 1000 * (nansum(W.totmf1neg.(reg)) ./ num_in_region);
    meridnet.(reg)  = 1000 * (nansum(W.totmf1pos.(reg) + W.totmf1neg.(reg)) ./ num_in_region);
    
    %frac tot
    fractot.(reg) = round(W.fractotmf.(reg),2);
    
    
%     disp({'Eastward Flux', round(1000*W.totmf2pos.A / (W.totmf2pos.An + W.totmf2neg.An),2),round(1000*W.totmf2pos.B / (W.totmf2pos.Bn + W.totmf2neg.Bn),2),round(1000*(W.totmf2pos.A + W.totmf2neg.A + W.totmf2pos.B + W.totmf2neg.B) / (W.totmf2pos.An + W.totmf2neg.An + W.totmf2pos.Bn + W.totmf2neg.Bn),2)})
%     disp({'Westward Flux', round(1000*W.totmf2neg.A / (W.totmf2pos.An + W.totmf2neg.An),2),round(1000*W.totmf2neg.B / (W.totmf2pos.Bn + W.totmf2neg.Bn),2),round(1000*(W.totmf2pos.A + W.totmf2neg.A + W.totmf2pos.B + W.totmf2neg.B) / (W.totmf2pos.An + W.totmf2neg.An + W.totmf2pos.Bn + W.totmf2neg.Bn),2)})
%     disp({'Northward Flux',round(1000*W.totmf1pos.A / (W.totmf1pos.An + W.totmf1neg.An),2),round(1000*W.totmf1pos.B / (W.totmf1pos.Bn + W.totmf1neg.Bn),2),round(1000*(W.totmf1pos.A + W.totmf1neg.A + W.totmf1pos.B + W.totmf1neg.B) / (W.totmf1pos.An + W.totmf1neg.An + W.totmf1pos.Bn + W.totmf1neg.Bn),2)})
%     disp({'Southward Flux',round(1000*W.totmf1neg.A / (W.totmf1pos.An + W.totmf1neg.An),2),round(1000*W.totmf1neg.B / (W.totmf1pos.Bn + W.totmf1neg.Bn),2),round(1000*(W.totmf1pos.A + W.totmf1neg.A + W.totmf1pos.B + W.totmf1neg.B) / (W.totmf1pos.An + W.totmf1neg.An + W.totmf1pos.Bn + W.totmf1neg.Bn),2)})
%     disp('     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%     disp({'Mean Flux',round(1000*nansum(W.totmf.A)./nansum(W.totmf.An),2),round(1000*nansum(W.totmf.B)./nansum(W.totmf.Bn),2),round(1000*nansum(W.totmf.C)./nansum(W.totmf.Cn),2)})
%     disp({'Total Flux %',round(W.fractotmf.A,2),round(W.fractotmf.B,2),[]})
%     disp({'Gini Coeff',round(W.gini.A,2),round(W.gini.B,2),round(W.gini.C,2)})
%     
    end
    
    % DISPLAY!!
    
    disp(['%%%%%%%%%%% ' names{w} ' %%%%%%%%%%%%%%%%%%%%%%%%%'])
    
    disp({'Region           ',' A ',' B ',' C '})
    disp(' ')
    disp({'Mean Amplitude   ',round(meanamp.A,2),round(meanamp.B,2),round(meanamp.C,2)})
    disp(' ')
    disp({'Eastward Flux:   ',round(eastflux.A,2),round(eastflux.B,2),round(eastflux.C,2)})
    disp({'Westward Flux:   ',round(westflux.A,2),round(westflux.B,2),round(westflux.C,2)})
    disp({'Net Zonal Flux:  ',round(zonalnet.A,2),round(zonalnet.B,2),round(zonalnet.C,2)})
    disp(' ')
    disp({'Northward Flux:  ',round(northflux.A,2),round(northflux.B,2),round(northflux.C,2)})
    disp({'Southward Flux:  ',round(southflux.A,2),round(southflux.B,2),round(southflux.C,2)})
    disp({'Net Merid Flux:  ',round(meridnet.A,2),round(meridnet.B,2),round(meridnet.C,2)})
    disp(' ')
    disp({'Frac of Total MF:',round(fractot.A,2),round(fractot.B,2),round(fractot.C,2)})
    
    
    
    
end




% 
% disp('%%%%%%%%%%% MODEL %%%%%%%%%%%%%%%%%%%%%%%%')
% disp({'Mean Amp',round(nansum(M.a.A)./nansum(M.a.An),2),round(nansum(M.a.B)./nansum(M.a.Bn),2),round(nansum(M.a.C)./nansum(M.a.Cn),2)})
% disp({'Mean Flux',round(1000*nansum(M.totmf.A)./nansum(M.totmf.An),2),round(1000*nansum(M.totmf.B)./nansum(M.totmf.Bn),2),round(1000*nansum(M.totmf.C)./nansum(M.totmf.Cn),2)})
% disp({'Total Flux %',round(M.fractotmf.A,2),round(M.fractotmf.B,2)})
% disp({'Gini Coeff',round(M.gini.A,2),round(M.gini.B,2),round(M.gini.C,2)})
% 
% disp('%%%%%%%%%%% MODEL AS AIRS %%%%%%%%%%%%%%%%')
% disp({'Mean Amp',round(nansum(MA.a.A)./nansum(MA.a.An),2),round(nansum(MA.a.B)./nansum(MA.a.Bn),2),round(nansum(MA.a.C)./nansum(MA.a.Cn),2)})
% disp({'Mean Flux',round(1000*nansum(MA.totmf.A)./nansum(MA.totmf.An),2),round(1000*nansum(MA.totmf.B)./nansum(MA.totmf.Bn),2),round(1000*nansum(MA.totmf.C)./nansum(MA.totmf.Cn),2)})
% disp({'Total Flux %',round(MA.fractotmf.A,2),round(MA.fractotmf.B,2)})
% disp({'Gini Coeff',round(MA.gini.A,2),round(MA.gini.B,2),round(MA.gini.C,2)})
% 
% 





return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_Regions','A','M','MA','Regions','L','h','zlims')
disp('Saved.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE 4DP?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('/Users/neil/Drive/MATLAB/SGWEX/AirsModel_Regions_4dp','A','M','MA','Regions','L','h','zlims')
disp('Saved.')

return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTY WOTTY PLOT PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.525 0.6])

%-------------------------------------------------------
vert_gap    = 0.05;    horz_gap    = 0.05;
lower_marg  = 0.05;    upper_marg  = 0.05;
left_marg   = 0.05;     right_marg  = 0.05;

rows = 1; cols = 1;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 20;

%--------------------------------------------------------

load('/Users/neil/Drive/MATLAB/SGWEX/20150705-1700_Model_as_AIRS_3DST.mat')
load('/Users/neil/Drive/MATLAB/SGWEX/20150705-1700_AirsSG_3DST.mat')

axx = gca;

% % some AIRS perturbations?
% zlev = 27;
% tbp = (Model.Tpg_sc(:,:,zlev) + Airs.Tpg_sc(:,:,zlev)) ./ 2;
% % tbp = -(tbp);
% % hold on; pcolor(D2(:,:,zlev),D1(:,:,zlev),Airs.Tpg_sc(:,:,zlev)); shat;
% % hold on; contourf(D2(:,:,zlev),D1(:,:,zlev),Airs.Tpg_sc(:,:,zlev),-5:0.5:5,'edgecolor','none'); shat;
% hold on; imagesc(linspace(-600,600,80),linspace(-450,450,60),tbp); ydir;
% 
% clim([-5 5])
% colormap(gca,cbrew)

ynudge = 0;

% Map and topography:
hMap = surface(xvec,yvec+ynudge,tpr,Mapr,...
    'FaceColor','texturemap',...
    'EdgeColor','none',...
    'CDataMapping','direct');

% the ring...
theta = reverse(180:360);
x = r .* sind(theta) + (-L/2 + d1 + r);
y = r .* cosd(theta);

% outline
xco = [(-L/2) (L/2) (L/2) (-L/2) (-L/2)];
yco = [(h/2) (h/2) (-h/2) (-h/2) (h/2)];

% region A patch
xca = [(-L/2) (-L/2 + r + d1) x (-L/2 + r + d1) (-L/2) (-L/2)];
yca = [(h/2) (h/2) y (-h/2) (-h/2) (h/2)];

% region B patch
xcb = [(L/2) x (L/2) (L/2)];
ycb = [(h/2) y (-h/2) (h/2)];

% 250km ring...
theta250 = 0:360;
x250 = 250 .* sind(theta250) + 100;
y250 = 250 .* cosd(theta250);

%%%% OVERLAYS
overlays = 0;

if overlays

% Region A
hold on; P.A = patch(xca,yca,10.*ones(size(xca)),mcolor(6));
P.A.FaceAlpha = 0.35;
P.A.LineStyle = 'none';

% Region B
hold on; P.B = patch(xcb,ycb,10.*ones(size(xcb)),[1 0.6 0]);
P.B.FaceAlpha = 0.35;
P.B.LineStyle = 'none';

% Outlines...
hold on; plot3(xco,yco,11.*ones(size(xco)),'linest','--','color','w','linewi',3);
hold on; plot3(x,y,11.*ones(size(x)),'linest','--','color','w','linewi',3);
% hold on; plot(x250,y250,'linest','--','color',mcolor(7),'linewi',3); grid on;

% BIG letters:
hold on; nph_text([0.2 0.775],'A','fontsize',2*fs,'color','w');
hold on; nph_text([0.7 0.775],'B','fontsize',2*fs,'color','w');

end

% South Georgia coastline
clat = -54.5269874182436;
clon = -37.1367495437688;
LonBox = [-40 -30];
LatBox = [-60 -40];
C = nph_draw_coastline([LonBox(1) LatBox(1) ; LonBox(2) LatBox(2)],0,0,'noplot','color','k');
for i = 1:length(C)
    if length(C(i).Lon) > 1
        [d,az] = distance(clat,clon,C(i).Lat,C(i).Lon);
        hold on; plot3(deg2km(d).*sind(az),deg2km(d).*cosd(az)+ynudge,ones(size(d)),'color','y','linewi',2.5);
        hold on; plot3(deg2km(d).*sind(az),deg2km(d).*cosd(az)+ynudge,ones(size(d)),'color','k','linewi',1.5);
    end
end

% limits and ticks...
xlim([-600 600])
ylim([-400 400])

set(gca,'linewi',1.5,'layer','top','tickdir','out');

axx.XMinorTick = 'on';
axx.XAxis.MinorTickValues = -600:100:600;

xlabel('x (km)')
ylabel('y (km)')


% grid on


% Extra axes for lat lon ticks?
axx2 = axes('position',axx.Position,'color','none','yaxislocation','right','xaxislocation','top');
axx2.YLim = minmax(LAT);
axx2.XLim = minmax(LON);
axx2.LineWidth = 1.5;
axx2.Layer = 'top';
axx2.TickDir = 'out';
axx2.XTick = -60:5:-20;
axx2.YTick = -60:2:-30;
axx2.XTickLabel = lonlabels(axx2.XTick);
axx2.YTickLabel = latlabels(axx2.YTick);
% axx2.XMinorTick = 'on';
% axx2.XAxis.MinorTickValues = -60:-20;

setfont(fs)


return





%% EXPORT? ================================================================

savename = ['~/Desktop/AirsModel_Regions2'];

disp(['Exporting to ' savename '...'])

nph_saveas(gcf,savename,'png')

disp('Done.')







return






