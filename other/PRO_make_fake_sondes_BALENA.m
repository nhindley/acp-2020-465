
%
% designed for use on BALENA
%
% version 2.1 - 4D gridded interpolant between model timesteps
%
% EDIT: NEW VERSION to use linear interpolation between timesteps
%
% Just about fits in memory, but going to have to use BALENA due to
% accessing the model timesteps.


% %% PARALLELISATION SETTINGS:
% % delete(gcp('nocreate'));
% 
% % maxNumCompThreads(16);
% % parpool(16);


%% LOAD RADIOSONDE FILES ==================================================
% type = '2'; %sec

matlabdirec = '/home/f/nh351/MATLAB/';
modeldirec = '/home/f/nh351/scratch/sgwex/';
sondedirec = '/home/f/nh351/scratch/sgwex/sondes/';
savedirec = '/home/f/nh351/scratch/sgwex/fakesondes/';

sondefiles = dir([sondedirec '*Sonde*.mat']);

for i = length(sondefiles):-1:1
    
    Mu = []; Mv = []; Mw = []; Mt = [];
    
    %% SONDE ===================================================================
    tic
    disp(['%%% ' num2str(i) ' of ' num2str(length(sondefiles)) ' %%% -  Sonde ' sondefiles(i).name ' ...'])
    
    S = load([sondedirec sondefiles(i).name]);
    
    % project sonde wind into u and v:
    S.Sonde.u = S.Sonde.WindSpeed.*sind(S.Sonde.WindDirection);
    S.Sonde.v = S.Sonde.WindSpeed.*cosd(S.Sonde.WindDirection);
    
    % track Sonde's X and Y positions
    [arclen,az] = distance(S.Sonde.Lat(1),S.Sonde.Lon(1),S.Sonde.Lat,S.Sonde.Lon);
    d = deg2km(arclen);
    
    S.Sonde.x = d.*sind(az); % ZONAL
    S.Sonde.y = d.*cosd(az); % MERIDIONAL
    S.Sonde.z = S.Sonde.Alt;
    
    %% MODEL ===============================================================
    
    
    % Sonde flight will be spread across a few Model time steps.
    % find how many model timesteps this sonde traversed:
    desired_timesteps = unique(floor(S.Sonde.Time*24))./24;
    
    % also include the extra one at the end to go to:
    desired_timesteps(end+1) = (desired_timesteps(end) + 1/24);
    
    % load all those model timesteps into a gridded interpolant:
    missing_timesteps = ones(size(desired_timesteps));
    for n = 1:length(desired_timesteps)
        
        timestep = datestr(desired_timesteps(n),'yyyymmdd-HHMM');
        
        % find the runname for this timestep:
        runname = 'not_assigned'; formatin = 'yyyymmdd-HHMM';
        if inrange(datenum(timestep,formatin),datenum('20150702-1300',formatin),datenum('20150707-1200',formatin))
            runname = 'u-ab326';
        end
        if inrange(datenum(timestep,formatin),datenum('20150612-1300',formatin),datenum('20150702-1200',formatin))
            runname = 'u-ae766';
        end
        if inrange(datenum(timestep,formatin),datenum('20150101-1300',formatin),datenum('20150131-1200',formatin))
            runname = 'u-ag477';
            if inrange(datenum(timestep,formatin),datenum('20150108-1300',formatin),datenum('20150109-1200',formatin))
                runname = 'not_assigned'; % avoid the missing data on day 20150108.
            end
        end
        if inrange(datenum(timestep,formatin),datenum('20130701-1300',formatin),datenum('20130731-2300',formatin))
            runname = 'dlmoa';
        end
        
        % loady load...
        try % in case a sonde was launched outside of the model window:
        
            disp(['Loading Model timestep ' timestep '...'])
            Model = load([modeldirec runname '/' timestep '_SG_uvwtpd.mat']);
            
            % Store all model timesteps in 4D object
            Mu = cat(4,Mu,Model.Model.u);
            Mv = cat(4,Mv,Model.Model.v(:,1:end-1,:)); % fix offset v grid
            Mw = cat(4,Mw,Model.Model.w);
            Mt = cat(4,Mt,Model.Model.temp);
        
        catch err
            disp(['COULDN''T LOAD MODEL TIMESTEP ' timestep ' from run ''' runname '''.'])
            missing_timesteps(n) = 0;
        end
        
    end
    
    desired_timesteps = desired_timesteps(logical(missing_timesteps));
    
    % If we couldn't find any model timesteps for the launch, ignore this
    % sonde.
    if length(desired_timesteps) <= 1
        continue
    end
    
    
    %     % fix offset grid
    %     Model.v = Model.v(:,1:end-1,:);
    
    % % ndgrid required for griddedInterpolant
    %     [D1,D2,D3,D4] = ndgrid(...
    %         single(linspace(-450,450,600)),...
    %         single(linspace(-600,600,800)),...
    %         single(Model.Alt./1000),...
    %         single(desired_timesteps)); % MODEL ALT IS IN METRES!!!!
    %
    
    disp('Creating gridded interpolant...')
    
    % remember start timestep for singling the data - matlab will loose
    % precision and get non-unique grid vectors. Which apparently is a problem.
    start_timestep = min(desired_timesteps);
    
    % use grid vectors rather than ndgrid or meshgrid. Saves on time and
    % memory. (Wow seems to make a big difference!)
    grid_vectors = {single(linspace(-600,600,800)),...
        single(linspace(-450,450,600)),...
        single(Model.Model.Alt./1000),...
        single(desired_timesteps-start_timestep)};
    
    % interolants:
    Fu = griddedInterpolant(grid_vectors,Mu,'linear','none');
    Fv = griddedInterpolant(grid_vectors,Mv,'linear','none');
    Fw = griddedInterpolant(grid_vectors,Mw,'linear','none');
    Ft = griddedInterpolant(grid_vectors,Mt,'linear','none');
    
    % evaluate:      % note the change in x and y
    FSu = Fu(S.Sonde.x,S.Sonde.y,S.Sonde.z,S.Sonde.Time-start_timestep);
    FSv = Fv(S.Sonde.x,S.Sonde.y,S.Sonde.z,S.Sonde.Time-start_timestep);
    FSw = Fw(S.Sonde.x,S.Sonde.y,S.Sonde.z,S.Sonde.Time-start_timestep);
    FSt = Ft(S.Sonde.x,S.Sonde.y,S.Sonde.z,S.Sonde.Time-start_timestep);
    
    
    %% FAKE SONDES =============================================================
    
    % get type:
    type = sprintf('%02d',round(diff(S.Sonde.Time(1:2))*24*60*60));
    
    FakeSonde = struct;
    FakeSonde.StartTime = S.Sonde.StartTime;
    
    FakeSonde.Meta = struct;
    FakeSonde.Meta.Type = [type 's'];
    FakeSonde.Meta.CreatedOn = now;
    FakeSonde.Meta.CreatedBy = [mfilename '.m'];
    
    FakeSonde.Time      = S.Sonde.Time;
    FakeSonde.Lat       = S.Sonde.Lat;
    FakeSonde.Lon       = S.Sonde.Lon;
    FakeSonde.Alt       = S.Sonde.Alt;
    FakeSonde.u         = FSu;
    FakeSonde.v         = FSv;
    FakeSonde.w         = FSw;
    FakeSonde.Temp      = FSt;
    
    FakeSonde.Meta.Units = {
        'Time (Matlab Time, decimal days since 00-JAN-0000)'
        'Altitude (km)'
        'Temperature (K)'
        'Wind Speed (m/s)'};
    
    savename = [sondefiles(i).name(1:13) '_' type 's_FakeSonde.mat'];
    
%     parsave([savedirec savename],FakeSonde,'FakeSonde');
    save([savedirec savename],'FakeSonde');
    
    tim = toc;
    disp(['Saved ' savename ' in ' num2str(round(tim)) 's.'])
    
    
    
end









return











%% LOAD RADIOSONDE FILES ==================================================
type = '10'; %sec

direc = ['/Users/neil/data/SGWEX-Radiosondes/' type 'second/matlab/'];

% WINTER CAMPAIGN FOR NOW

campaign = '2015';

files = cat(1,dir([direc campaign '06*.mat']),dir([direc campaign '07*.mat']));


d = [];

for i = 1:length(files)
    
    disp(['%%% ' num2str(i) ' of ' num2str(length(files)) ' %%% -  Loading Sonde ' files(i).name ' ...'])
    
    load([direc files(i).name]);
    
    % get timestep
    % round to next hour (fair for time of flight)
    granule = datestr(dateshift(datetime(S.Sonde.StartTime),'end','hour'),'yyyymmdd-HHMM');
    
    % project sonde wind into u and v:
    S.Sonde.u = S.Sonde.WindSpeed.*sind(S.Sonde.WindDirection);
    S.Sonde.v = S.Sonde.WindSpeed.*cosd(S.Sonde.WindDirection);
    
    % track Sonde's X and Y positions
    [arclen,az] = distance(S.Sonde.Lat(1),S.Sonde.Lon(1),S.Sonde.Lat,S.Sonde.Lon);
    d = deg2km(arclen);
    
    S.Sonde.x = d.*sind(az);
    S.Sonde.y = d.*cosd(az);
    S.Sonde.z = S.Sonde.Alt;
    
    %% LOAD MODEL =============================================================
    
    disp(['Loading Model ' granule ' ...'])
    
    % DUMMY MODEL: Hard drive is busy!!
    %     load('/Users/neil/data/SGWEX/model3dst/SG/20150705-1700_SG_uvwtpd.mat')
    try
        switch campaign
            case '2015'
                load(['/Volumes/RED/Model/SG/' granule '_SG_uvwtpd.mat'])
            case '2013'
                load(['/Volumes/BLUE/Model/SG/' granule '_SG_uvwtpd.mat'])
        end
        
        
        % fix offset grid
        Model.v = Model.v(:,1:end-1,:);
        
        % ndgrid required for griddedInterpolant
        [D1,D2,D3] = meshgrid(...
            single(linspace(-450,450,600)),...
            single(linspace(-600,600,800)),...
            single(Model.Alt./1000)); % MODEL ALT IS IN METRES!!!!
        
        D1 = permute(D1,[2 1 3]);
        D2 = permute(D2,[2 1 3]);
        D3 = permute(D3,[2 1 3]);
        
        disp('Evaluating sonde path through gridded model interpolant...')
        
        Fu = griddedInterpolant(D1,D2,D3,permute(Model.u,[2 1 3]),'nearest','none');
        Modelu = Fu(S.Sonde.y,S.Sonde.x,S.Sonde.z);
        
        Fv = griddedInterpolant(D1,D2,D3,permute(Model.v,[2 1 3]),'nearest','none');
        Modelv = Fv(S.Sonde.y,S.Sonde.x,S.Sonde.z);
        
        Ft = griddedInterpolant(D1,D2,D3,permute(Model.temp,[2 1 3]),'nearest','none');
        Modelt = Ft(S.Sonde.y,S.Sonde.x,S.Sonde.z);
        
        Fp = griddedInterpolant(D1,D2,D3,permute(Model.pres,[2 1 3]),'nearest','none');
        Modelp = Fp(S.Sonde.y,S.Sonde.x,S.Sonde.z);
        
        %% PLOT ===============================================================
        
        figure;
        hold all; whitefig; grid on;
        
        hold on; plot(S.Sonde.u,S.Sonde.Alt,'r');
        hold on; plot(Modelu,S.Sonde.Alt,'b');
        
        hold on; plot(S.Sonde.v,S.Sonde.Alt,'r');
        hold on; plot(Modelv,S.Sonde.Alt,'b');
        
        drawnow;
        return
        %% SAVE ===============================================================
        
        FakeSonde = struct;
        %     if strcmpi(sondefiles())
        FakeS.Sonde.Type = [type 's'];
        FakeS.Sonde.StartTime = S.Sonde.StartTime;
        FakeS.Sonde.ModelTimestep = granule;
        
        FakeS.Sonde.Time      = S.Sonde.Time;
        FakeS.Sonde.Alt       = S.Sonde.Alt;
        FakeS.Sonde.u         = Modelu;
        FakeS.Sonde.v         = Modelv;
        FakeS.Sonde.Temp      = Modelt;
        FakeS.Sonde.Pressure  = Modelp;
        
        FakeS.Sonde.Lat = S.Sonde.Lat;
        FakeS.Sonde.Lon = S.Sonde.Lon;
        
        FakeS.Sonde.Units = {
            'Time (Matlab Time, decimal days since 00-JAN-0000)'
            'Altitude (km)'
            'Pressure (hPa)'
            'Temperature (K)'
            'Wind Speed (m/s)'};
        
        savedirec = ['/Users/neil/data/SGWEX-Radiosondes/' type 'second/matlab/FakeSondes/'];
        save([savedirec files(i).name(1:13) '_' type 's_FakeS.Sonde.mat'],'FakeSonde');
        
        disp(['Saved ' files(i).name(1:13) '_' type 's_FakeS.Sonde.mat'])
        
        
    catch err
        disp(err.message)
    end
    
    
end



return


%%




%% EXPERIMENT WITH 4-D GRIDDED INTERPOLANT:

sondehours = unique(floor(S.Sonde.Time .* 24) ./ 24,'stable');
sondehours(end+1) = sondehours(end) + 1/24; % also need the next hour to interpolate to

[D1,D2,D3,D4] = ndgrid(...
    single(linspace(-450,450,600)),...
    single(linspace(-600,600,800)),...
    single(Model.Alt./1000),...
    single(1:4));


%%

m = permute(Model.u,[2 1 3]);

ur = cat(4,rand+m,rand+m,rand+m,rand+m);

%%

% ur = repmat(,1,1,1,length(sondehours));

Fu = griddedInterpolant(D1,D2,D3,D4,ur,'nearest','none');

Modelur = Fu(S.Sonde.y,S.Sonde.x,S.Sonde.z,S.Sonde.Time);


%%
figure; hold all; whitefig;

hold on; plot(S.Sonde.u,S.Sonde.Alt,'r'); grid on;
hold on; plot(Modelur,S.Sonde.Alt,'b'); grid on;


%%
figure; hold all; whitefig;

hold on; plot(Modelu,S.Sonde.Alt,'r'); grid on;
hold on; plot(Modelur,S.Sonde.Alt,'b'); grid on;




































