

% Process SGWEX South Georgia Radiosonde data into matlab format.

type = '2';

direc = ['/Users/neil/data/SGWEX-Radiosondes/' type 'second/'];
% direc = '/Volumes/NJMitchell-Scratch/Data/SGWEX-Radiosondes/2second/';
% direc = '/Users/neil/data/SGWEX-Radiosondes/10second/';

files = dir([direc '*.dat']);

% for i = 108:108
for i = 1:length(files)
    
    try
        
        Sonde = getsonde([direc files(i).name],type);
        
        savename = [datestr(datenum(Sonde.StartTime),'yyyymmdd-HHMM') '_' type 's_Sonde.mat'];
        save([direc 'matlab/' savename],'Sonde');
        
        disp([files(i).name ' converted and saved!'])
        
    catch err
        
        disp([files(i).name ' failed: ' err.message])

    end


end








