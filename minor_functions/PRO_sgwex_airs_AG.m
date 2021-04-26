
%%
% this is a launcher for the model-as-airs_AG processor...

% Generate and execute a slurm script to run the MODEL AS AIRS processing
% for a given day. Doing day by day in separate job each should be much
% faster overall :)


startup

% % day range for both 2013 and 2015 model runs
% dayrange = [...
%     datenum(2013,07,01):datenum(2013,08,01) ...
%     datenum(2015,06,12):datenum(2015,07,06) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DIRECTORIES

% for loading...
airsdirec = '/home/f/nh351/scratch/data/sgwex/Airs/three_corner_overlapping_granule_pairs/matlab/';

% and for saving...
savedirec = '/home/f/nh351/scratch/data/sgwex/3DST/Airs/AG_sm/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each day, submit a new job...
% for d = 1:length(dayrange)
% for d = 56
      
%     yyyymmdd = datestr(dayrange(d),'yyyymmdd');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% GENERATE SLURM SCRIPT ===============================================    
    Slurm = struct;
    
    Slurm.Name = 'AirsSG_3DST_AG';
%     Slurm.Batch = 'batch-short';
    Slurm.Batch = 'batch-all';
    
    % choose max time allocation for balena job:
    % default 360 minutes. Might get nudged up the queue if you reduce this.
    Slurm.WallTime = 60;
    
    Slurm.Filename = ['/home/f/nh351/MATLAB/slurm/scripts/' Slurm.Name '.slurm'];
    
    Slurm.Command = '';
    
    %%%% SLURM COMMANDS
    Slurm.Command = [Slurm.Command ...
        'airsdirec = ''' airsdirec '''; ' ...
        'savedirec = ''' savedirec '''; ' ...
        'FUN_sgwex_airs_regrid_and_3dst_AG(airsdirec,savedirec)'];
    
    % NOTE: If you get a problem with parpool() "not enough input
    % arguments", on balena go to ~/.matlab and delete the folder called
    % local_cluster_jobs. This will probably terminate any current jobs
    % too, so check first.
    
    
    % finally, write the slurm script:
    writeslurm(Slurm);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% SUBMIT AS BALENA JOB
    
    disp('==============================================================')
    disp(['Submitting slurm script ' Slurm.Filename '...'])
    
    system(['sbatch ' Slurm.Filename]);
    
    
    
% end




return




% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % % Generate and execute a slurm script to 3DST and then regrid the 3DST
% % % % % output, or do one or the other depending on what you want to do.
% % % % 
% % % % % NEW PLAN, set up and queue a slurm script for each AIRS 3D day. You can
% % % % % get about 7 days done in one six hour job, but it's probably better to to
% % % % % just submit lots of one-day jobs for now, particularly to batch-short,
% % % % % which seems to be pretty efficient when used like this.
% % % % 
% % % % startup
% % % % 
% % % % years = 2009:2010;
% % % % 
% % % % % year:
% % % % for yr = years
% % % %     
% % % %     % yr = 2003;
% % % %     
% % % %     % dayrange:
% % % %     dayrange = daynumber(datenum(yr,01,01):datenum(yr,12,31));
% % % %     % dayrange = daynumber(datenum(yr,01,01));
% % % %     
% % % %     
% % % %     % =========================================================================
% % % %     % What do you want to do? 3DST and regrid or regrid only or...?
% % % %     st = 0;
% % % %     regrid = 1;
% % % %     useparpool = 0;
% % % %     
% % % %     if st == 1 && regrid == 0
% % % %         name = '3dst';
% % % %     end
% % % %     if regrid == 1 && st == 0
% % % %         name = 'regrid';
% % % %     end
% % % %     if st == 1 && regrid == 1
% % % %         name = '3dst_and_regrid';
% % % %     end
% % % %     
% % % %     % =========================================================================
% % % %     
% % % %     
% % % %     for d = 1:length(dayrange)
% % % %         
% % % %         % GENERATE SLURM SCRIPT ===============================================
% % % %         yyyymmdd = datestr(datenum(yr,00,dayrange(d)),'yyyymmdd');
% % % %         
% % % %         Slurm = struct;
% % % %         
% % % %         Slurm.Name = [yyyymmdd '_airs_' name '_z40km'];
% % % %         Slurm.Batch = 'batch-short';
% % % %         %         Slurm.Batch = 'batch-all';
% % % %         
% % % %         % choose max time allocation for balena job:
% % % %         WALLTIME = 10; % default 360 minutes. Might get nudged up the queue if you reduce this.
% % % %         Slurm.Walltime = WALLTIME;
% % % %         
% % % %         Slurm.Filename = ['/home/f/nh351/MATLAB/slurm/scripts/' Slurm.Name '.slurm'];
% % % %         
% % % %         Slurm.Command = '';
% % % %         
% % % %         Slurm.WallTime = WALLTIME; % wall time in minutes.
% % % %         
% % % %         % ---------------------------------------------------------------------
% % % %         % Write parallelisation and 3DST commands...
% % % %         if st == 1
% % % %             if useparpool
% % % %                 Slurm.Command = [Slurm.Command ...
% % % %                     'delete(gcp(''nocreate'')); maxNumCompThreads(16); parpool(''local'',16); ' ...
% % % %                     'yr = ' num2str(yr) '; dayrange = ' num2str(dayrange(d)) '; ' ...
% % % %                     'FUN_airs_3dst_z40km_BALENA(yr,dayrange); '];
% % % %             else
% % % %                 Slurm.Command = [Slurm.Command ...
% % % %                     'delete(gcp(''nocreate'')); yr = ' num2str(yr) '; dayrange = ' num2str(dayrange(d)) '; ' ...
% % % %                     'FUN_airs_3dst_z40km_BALENA(yr,dayrange); '];
% % % %             end
% % % %         end
% % % %         % NOTE: If you get a problem with parpool() "not enough input
% % % %         % arguments", on balena go to ~/.matlab and delete the folder called
% % % %         % local_cluster_jobs. This will probably terminate any current jobs
% % % %         % too, so check first.
% % % %         
% % % %         % ---------------------------------------------------------------------
% % % %         % Add a section to regrid the ST variables onto a regular 1x1 degree
% % % %         % lat lon grid for plotting later:
% % % %         if regrid == 1
% % % %             if useparpool
% % % %                 Slurm.Command = [Slurm.Command ...
% % % %                     'delete(gcp(''nocreate'')); maxNumCompThreads(16); parpool(''local'',16); ' ...
% % % %                     'yr = ' num2str(yr) '; dayrange = ' num2str(dayrange(d)) '; ' ...
% % % %                     'FUN_airs3d_regridder_z40km_BALENA(yr,dayrange);'];
% % % %             else
% % % %                 Slurm.Command = [Slurm.Command ...
% % % %                     'delete(gcp(''nocreate''));' ...
% % % %                     'yr = ' num2str(yr) '; dayrange = ' num2str(dayrange(d)) '; ' ...
% % % %                     'FUN_airs3d_regridder_z40km_BALENA(yr,dayrange);'];
% % % %                 % ^ remember to kill the parpool first, won't need it.
% % % %             end
% % % %         end
% % % %         writeslurm(Slurm);
% % % %         
% % % %         %     disp(['Saved ' yyyymmdd '_airs_3dst_z40km.slurm'])
% % % %         
% % % %         
% % % %         
% % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %         %% SUBMIT AS BALENA JOB
% % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %         
% % % %         disp('==============================================================')
% % % %         disp(['Submitting slurm script ' Slurm.Filename '...'])
% % % %         
% % % %         system(['sbatch ' Slurm.Filename]);
% % % %         
% % % %         
% % % %         
% % % %         
% % % %         
% % % %     end % next day
% % % %     
% % % % end % next year
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % %
% % % % % %% PARALLELISATION SETTINGS:
% % % % % % delete(gcp('nocreate')) % delete any currently active session.
% % % % % maxNumCompThreads(16);
% % % % % parpool(16);
% % % % %
% % % % % % Load and process a load of 3D Airs with the 3DST to see if we can
% % % % % % reproduce the bias Corwin has found. If not, well, that just leaves more
% % % % % % questions...
% % % % %
% % % % % yr = 2015;
% % % % % yrstr = num2str(yr);
% % % % %
% % % % % dayrange = 213:222;
% % % % % dayrange = 223:233;
% % % % % dayrange = 233:243;
% % % % %
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % on_balena = 1;
% % % % % % DIRECS
% % % % % switch on_balena
% % % % %     case 1
% % % % %         datadirec = '/home/f/cw785/scratch/Data/AIRS/airs3d/';
% % % % %         savedirec = '/home/f/nh351/scratch/data/airs3d/3DST/z40km/';
% % % % %     otherwise
% % % % %         datadirec = '/Volumes/SDCard/data/airs3d/';
% % % % %         savedirec = '/Volumes/SDCard/data/airs3d/3DST/z40km/';
% % % % % end
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %
% % % % % % FOR EACH DAY
% % % % % for d = 1:length(dayrange)
% % % % %
% % % % %     direc = [datadirec  yrstr '/' sprintf('%03d',dayrange(d)) '/'];
% % % % %
% % % % %     files = dir([direc 'airs*.nc']);
% % % % %
% % % % %     % Storage:
% % % % %     A = nan(90,135,240);
% % % % %     BG = nan(90,135,240);
% % % % %     K = nan(90,135,240);
% % % % %     L = nan(90,135,240);
% % % % %     M = nan(90,135,240);
% % % % %     LAT = nan(90,135,240);
% % % % %     LON = nan(90,135,240);
% % % % %
% % % % %     % FOR EACH GRANULE
% % % % %     parfor i = 1:length(files)
% % % % %
% % % % %         disp([num2str(i) ' of ' num2str(length(files)) ' ----- ' files(i).name])
% % % % %
% % % % %         try
% % % % %
% % % % %             tic
% % % % %
% % % % %             % load:
% % % % %             nc1 = getnet([direc files(i).name]);
% % % % %
% % % % %             %% Extract variables:
% % % % %             T = nc1.Data.ret_temp;
% % % % %             T = permute(T,[2 3 1]);
% % % % %             Lat = nc1.Data.l1_lat;
% % % % %             Lon = nc1.Data.l1_lon;
% % % % %             Alt = nc1.Data.ret_z;
% % % % %
% % % % %             % NaN checker:
% % % % %             T(isnan(T)) = nanmean(T(:));
% % % % %
% % % % %             % remove bg:
% % % % %             [Tp,bg] = sgwex_airs_remove_bg(T);
% % % % %
% % % % %             % regular horz grid:
% % % % %             [Long,Latg,Tpg,xt_spacing,at_spacing] = put_airs_on_regular_grid(Lon,Lat,Tp);
% % % % %             [~,~,bgg,~,~] = put_airs_on_regular_grid(Lon,Lat,bg);
% % % % %
% % % % %             % regular alt grid:
% % % % %             vert_spacing = 4;
% % % % %             desired_alt_range = 8:vert_spacing:68;
% % % % %             Altg = desired_alt_range;
% % % % %
% % % % %             FT = griddedInterpolant({1:90,1:135,Alt},Tpg,'linear','none');
% % % % %             Tpg = FT({1:90,1:135,desired_alt_range});
% % % % %
% % % % %             Fbg = griddedInterpolant({1:90,1:135,Alt},bgg,'linear','none');
% % % % %             bgg = Fbg({1:90,1:135,desired_alt_range});
% % % % %
% % % % %             % remove amplitude increase with height:
% % % % %             ref_z = 40; % km
% % % % %             H = scale_height(ref_z);
% % % % %             sfvec = (exp(-(Altg-ref_z) ./ (2*H)))';
% % % % %             sf = permute(repmat(sfvec,1,90,135),[2 3 1]);
% % % % %             Tpg_sc = Tpg .* sf;
% % % % %
% % % % %             % Apply altitude window to roughly 20-60km:
% % % % %             winvec = [0 0 tukeywin(13,0.4)' 0]';
% % % % %             win = permute(repmat(winvec,1,90,135),[2 3 1]);
% % % % %             Tpg_sc = Tpg_sc .* win;
% % % % %
% % % % %             % apply smoother
% % % % %             Tpg_sc = smoothn(Tpg_sc,[3 3 1]);
% % % % %             bgg = smoothn(bgg,[3 3 1]);
% % % % %
% % % % %             t0 = toc;
% % % % %
% % % % %             disp(['Loaded and regridded in ' num2str(round(t0)) 's.'])
% % % % %
% % % % %             %% NDST ===============================================================
% % % % %
% % % % %             sampling_intervals = [xt_spacing at_spacing vert_spacing];
% % % % %             c = [0.25 0.25 0.25];
% % % % %
% % % % %             % calculate scales from selected wavelengths...
% % % % %             l1 = -1800:100:1800;
% % % % %             l2 = -2400:100:2400; % set for ONE granule ONLY.
% % % % %             l3 = -[60 30 20 15 12 10 8 6 5];
% % % % %
% % % % %             s1 = unique(round(xt_spacing .* size(Tpg_sc,1) ./ l1));   s1 = s1(isfinite(s1));
% % % % %             s2 = unique(round(at_spacing .* size(Tpg_sc,2) ./ l2));   s2 = s2(isfinite(s2));
% % % % %             s3 = unique(round(vert_spacing .* size(Tpg_sc,3) ./ l3)); s3 = s3(isfinite(s3));
% % % % %
% % % % %             % disp('3DST...')
% % % % %             tic
% % % % %             ST = nph_ndst(Tpg_sc,{s1,s2,s3},sampling_intervals,c);
% % % % %             t1 = toc;
% % % % %
% % % % %             disp(['3DST computed in ' num2str(round(t1)) 's.'])
% % % % %
% % % % %
% % % % %             %% COMPUTE ANGLES AND PROJECT:
% % % % %
% % % % %             ST.kh = quadadd(ST.F1,ST.F2);
% % % % %
% % % % %             % find angle CLOCKWISE from along track direction...
% % % % %             ang_at = atan2d(ST.F1,ST.F2);
% % % % %
% % % % %             % find azimuth of AT direction CLOCKWISE from north...
% % % % %             xt_mid = 45;
% % % % %             [~,az_at] = distance(Latg(xt_mid,1:end-1),Long(xt_mid,1:end-1),Latg(xt_mid,2:end),Long(xt_mid,2:end));
% % % % %             az_at(end+1) = az_at(end);
% % % % %             az_at = repmat(az_at,90,1);
% % % % %
% % % % %             % add angles... (since they both should be clockwise from north)
% % % % %             ang_north = wrapTo180(az_at + ang_at);
% % % % %
% % % % %             % finally, find k and l:
% % % % %             ST.k = ST.kh .* sind(ang_north);
% % % % %             ST.l = ST.kh .* cosd(ang_north);
% % % % %
% % % % %             %% COLLECT ALL TODAY'S GRANULES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %             zlevel = 9; % should be 40km for 8:4:68km vertical grid
% % % % %
% % % % %             A(:,:,i)  = ST.HA(:,:,9);
% % % % %             BG(:,:,i) = bgg(:,:,9);
% % % % %             K(:,:,i)  = ST.k(:,:,9);
% % % % %             L(:,:,i)  = ST.l(:,:,9);
% % % % %             M(:,:,i)  = ST.F3(:,:,9);
% % % % %
% % % % %             LAT(:,:,i) = Latg;
% % % % %             LON(:,:,i) = Long;
% % % % %
% % % % %         catch err
% % % % %
% % % % %             disp(['Can''t process granule ' files(i).name ': ' err.identifier])
% % % % %
% % % % %         end % end try for this granule
% % % % %     end % next granule (PARFOR)
% % % % %
% % % % %
% % % % %     %% SAVE %%%%%%%%%%%%%%
% % % % %
% % % % %     Airs = struct;
% % % % %
% % % % %     Airs.Lon = LON;
% % % % %     Airs.Lat = LAT;
% % % % %     Airs.z = 40;
% % % % %
% % % % %     Airs.A = A;
% % % % %     Airs.bg = BG;
% % % % %     Airs.k = K;
% % % % %     Airs.l = L;
% % % % %     Airs.m = M;
% % % % %
% % % % %     % Convert to singles:
% % % % %     f = fieldnames(Airs);
% % % % %     for i = 1:length(f)
% % % % %         Airs.(f) = single(Airs.(f));
% % % % %     end
% % % % %
% % % % %     savename = [datestr(datenum(yr,00,dayrange(d)),'yyyymmdd') '_AIRS_3DST_z40km.mat'];
% % % % %     %         parsave([savedirec savename],Day,'Day');
% % % % %     save([savedirec savename],'Day');
% % % % %
% % % % %     disp(['Saved to ' savename '.'])
% % % % %
% % % % %
% % % % % end % next day
% % % % %
% % % % %
% % % % %
% % % % % return
% % % % %
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
