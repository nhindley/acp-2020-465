

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB STARTUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('Running Matlab startup file...')

drivepath = '/Users/neil/Drive/MATLAB/';

cd(drivepath)

% Add Google Drive path:
paths = strsplit(genpath(drivepath),':');
addpath(paths{:});

% Revert to old Figure toolbar:
try %#ok
    if ~verLessThan('matlab','9.5')
        set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig));
        set(groot,'defaultAxesCreateFcn',  @(ax,~)set(ax.Toolbar,'Visible','off'));
    end
end

set(0,'defaulttextfontsize',14);
set(0,'defaultaxesfontsize',14);


clear;

