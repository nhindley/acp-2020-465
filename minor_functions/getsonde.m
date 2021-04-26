
% UPGRADE of nph_read_sonde.m
%
% Quick function for reading the data (and only the data) out of
% HiResRadioSonde files like the ones located in:
%
% /Volumes/NJMitchell-Scratch/Data/SGWEX-Radiosondes/2second/
%
% No idea if it'll work on other file types. Probably not. I don't think
% it'll work on even the 10 second files in the same directory.
%
% Use [type] to switch between 2 second format and 10 second format
%
% NOTE: 10s sondes seem to have wind speed in KNOTS.

function Sonde = getsonde(filename,type)

if ~exist('type','var')
    type = '2';
else
    if isnumeric(type)
        type = num2str(type);
    end
end

switch type
    % 2 SECOND DATA
    case {'2 second','2s','2sec','2'}
        
        % Assume lauch at KEP
        StartLon = -36.5;
        StartLat = -54.283;
        
        %% READ FILE
        delimiter = ' ';
        formatSpec = '%q%q%q%q%q%q%q%q%q%q%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);
        fclose(fileID);
        
        
        %% CONCATENATE to string object
        s = mat2cell(char(dataArray{1}), ones(length(dataArray{1}), 1));
        for i = 2:length(dataArray)
            s = cat(2,s,mat2cell(char(dataArray{i}), ones(length(dataArray{i}), 1)));
        end
        
        
        %% FIND WHERE THE START TIME IS:
        % using the "Started" heading as a start point:
        % add trailing whitespace...
        len = num2str(length(s{1,1}));
        sstar = sprintf(['% ' len 's'],'');
        sstar(1:7) = 'Started';
        t = strtrim(s(strcmpi(s(:,1),sstar),:));
        
        dd = sprintf('%02s',t{3});
        mmm = t{4};
        yyyy = ['20' sprintf('%02s',t{5})];
        HHMM = sprintf('%05s',t{6});
        
        
        %% FIND WHERE THE NUMERIC DATA STARTS
        % using the minute (min) heading as a start point:
        % add trailing whitespace...
        len = num2str(length(s{1,1}));
        mmin = sprintf(['% ' len 's'],'');
        mmin(1:3) = 'min';
        startrow = find(strcmpi(s(:,1),mmin))+1;
        
        % vars =
        % 1   2   3   4 5   6    7   8    9   10  11
        % min sec m/s m hPa degC %RH degC deg m/s empty
        
        
        %% COLLECT DATA:
        
        Sonde = struct;
        
        Sonde.StartTime = datestr(datenum([dd mmm yyyy '-' HHMM],'ddmmmyyyy-HH:MM'));
        
        % compute time:
        mins = (1/24/60) * str2double(s(startrow:end,1)); % decimal days
        secs = (1/24/60/60) * str2double(s(startrow:end,2));
        
        Sonde.Time            = datenum(Sonde.StartTime) + mins + secs; % matlab time
        Sonde.AscentRate      = str2double(s(startrow:end,3));
        Sonde.Alt             = str2double(s(startrow:end,4))./1000; % alt in km
        Sonde.Pressure        = str2double(s(startrow:end,5));
        Sonde.Temp            = str2double(s(startrow:end,6)) + 273.15; % K
        Sonde.RH              = str2double(s(startrow:end,7)); %
        Sonde.DewPoint        = str2double(s(startrow:end,8)); % deg C
        Sonde.WindDirection   = str2double(s(startrow:end,9)); % azimuth? (deg)
        Sonde.WindSpeed       = str2double(s(startrow:end,10)); % m/s
        
        % CHANGE WIND DIRECTION SO THAT ITS THE DIRECTION WIND IS GOING,
        % INSTEAD OF WHERE ITS COMING FROM:
        Sonde.WindDirection = wrapTo360(Sonde.WindDirection+180);
        
        
        Sonde.Units = {
            'Time (Matlab Time, decimal days since 00-JAN-0000)'
            'Ascent Rate (m/s)'
            'Altitude (km)'
            'Pressure (hPa)'
            'Temperature (K)'
            'Relative Humidity (%)'
            'Dew Point (degC)'
            'Wind Direction (deg from North)'
            'Wind Speed (m/s)'};
        
                %% FIND LAT LON POSITION
        % take sonde to be a leaf on the breeze, consider wind direction to
        % be a bearing from north. Take cumulative X and Y distances
        % covered every timestep
        
        samplingfreq = 10; % secs
        
        % need to tidy up NaNs here:
        Sonde.WindSpeed(isnan(Sonde.WindSpeed)) = nanmean(Sonde.WindSpeed);
        Sonde.WindDirection(isnan(Sonde.WindDirection)) = nanmean(Sonde.WindDirection);
        
        % now smooth over those added means...
        Sonde.WindSpeed = smooth(Sonde.WindSpeed,11);
        Sonde.WindDirection = smooth(Sonde.WindDirection,11);
        
        % find distance travelled every timestep, then sum for total
        arclen = km2deg((Sonde.WindSpeed*samplingfreq)./1000);
        
        % reckon
        [latout,lonout] = reckon(0,0,arclen,Sonde.WindDirection);
        
        Sonde.Lat = [StartLat ; StartLat+cumsum(latout(1:end-1))];
        Sonde.Lon = [StartLon ; StartLon+cumsum(lonout(1:end-1))];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        % 10 SECOND DATA
    case {'10 second','10 sec','10s','10'}
        
        % Assume launch at KEP
        StartLon = -36.5;
        StartLat = -54.283;
        
        fileID = fopen(filename,'r');
        formatSpec = '%8f%8f%8f%8f%8f%8f%8f%8f%f%[^\n\r]';
        dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);
        fclose(fileID);
        
        % convert to matrix
        d = nan(length(dataArray{1}),length(dataArray)-1);
        % can't use cell2mat because last column is always non-double for some reason:
        for i = 1:length(dataArray)-1
            d(:,i) = dataArray{i};
        end
        
        % Take top line for date finding:
        topline = strsplit(num2str(d(1,1:5)),' ');
        
        % Get matlab datenum:
        dd = sprintf('%02s',topline{3});
        mmm = monthname(str2double(topline{2}),'mmm');
        yyyy = ['20' sprintf('%02s',topline{1})];
        HH = sprintf('%02s',topline{4});
        MM = sprintf('%02s',topline{5});
        
        %% Create output variable - modified!!!
        % trim the top line:
        d = d(2:end,:);
        
        % NaNs instead of -999
        d(d == -999) = NaN;
        
        % compute time:
        mins = (1/24/60) * d(:,1); % decimal days
        secs = (1/24/60/60) * d(:,2);
        
        Sonde = struct;
        
        Sonde.StartTime       = datestr(datenum([yyyy mmm dd '-' HH MM],'yyyymmdd-HHMM'));
        Sonde.Time            = datenum(Sonde.StartTime) + mins + secs; % matlab time
        Sonde.Pressure        = d(:,3); % hPa
        Sonde.Alt             = d(:,4)./1000; % km
        Sonde.Temp            = d(:,5) + 273.15; % K
        Sonde.RH              = d(:,6); % relative humidity
        Sonde.DewPoint        = d(:,7); % dew point in deg C
        Sonde.WindDirection   = d(:,8); % Azimuth
        Sonde.WindSpeed       = d(:,9) * 0.514; % knots to m/s
        
        % CHANGE WIND DIRECTION SO THAT ITS THE DIRECTION WIND IS GOING,
        % INSTEAD OF WHERE ITS COMING FROM:
        Sonde.WindDirection = wrapTo360(Sonde.WindDirection+180);
        
        Sonde.Units = {
            'Time (Matlab Time, decimal days since 00-JAN-0000)'
            'Altitude (km)'
            'Pressure (hPa)'
            'Temperature (K)'
            'Relative Humidity (%)'
            'Dew Point (degC)'
            'Wind Direction (deg from North)'
            'Wind Speed (m/s)'};
        
        %% FIND LAT LON POSITION
        % take sonde to be a leaf on the breeze, consider wind direction to
        % be a bearing from north. Take cumulative X and Y distances
        % covered every timestep
        
        samplingfreq = 10; % secs
        
        % need to tidy up NaNs here:
        %try for loop based local mean instead, ie last known value
        Sonde.WindSpeed(isnan(Sonde.WindSpeed)) = nanmean(Sonde.WindSpeed);
        Sonde.WindDirection(isnan(Sonde.WindDirection)) = nanmean(Sonde.WindDirection);
        
        % now smooth over those added means...
        Sonde.WindSpeed = smooth(Sonde.WindSpeed,11);
        Sonde.WindDirection = smooth(Sonde.WindDirection,11);
        
        % find distance travelled every timestep, then sum for total
        arclen = km2deg((Sonde.WindSpeed*samplingfreq)./1000);
        
        % reckon
        [latout,lonout] = reckon(0,0,arclen,Sonde.WindDirection);
        
        Sonde.Lat = [StartLat ; StartLat+cumsum(latout(1:end-1))];
        Sonde.Lon = [StartLon ; StartLon+cumsum(lonout(1:end-1))];
        
        
end


%% FIXING STUFF ===========================================================
% wow radiosonde data is poor

% Altitude
% omg sometimes 10s sonde alts just flip sign for no reason..
Sonde.Alt = abs(Sonde.Alt);
Sonde.Alt(find(diff(Sonde.Alt) <= 0)+1) = NaN;

% Temperature (now outputting in K)
Sonde.Temp(Sonde.Temp < 0) = NaN;

% Pressure
Sonde.Pressure(Sonde.Pressure < 0) = NaN;
switch type
    case {'2 second','2s','2sec','2'}
        Sonde.Pressure = sgolayfilt(Sonde.Pressure,2,51);
    case {'10 second','10 sec','10s','10'}
        
        Sonde.Pressure = sgolayfilt(Sonde.Pressure,2,21);
end







% OUT.Pressure(OUT.Pressure < 0) = NaN;
% 
% switch type
%     case {'2 second','2s','2sec','2'}
%         
%         % find where change in pressure in 2 secs is more than 10 hPa
%         badinds = find(abs(diff(OUT.Pressure) > 10));
%         OUT.Pressure(badinds+1) = NaN;
%         
%         OUT.Pressure(OUT.Pressure < 5) = NaN; % very unlikely to reach 5 hPa
%         
%     case {'10 second','10 sec','10s','10'}
%         
%         % find where change in pressure in 10 secs is more than 50 hPa
%         badinds = find(abs(diff(OUT.Pressure) > 50));
%         OUT.Pressure(badinds+1) = NaN;
%         
%         OUT.Pressure(OUT.Pressure < 5) = NaN; % very unlikely to reach 5 hPa
%         
% end
















