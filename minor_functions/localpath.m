
%% DETERMINE WHAT PC WE ARE RUNNING ON IN ORDER TO GET PATHS CORRECT
% should have done this years ago

function path = localpath(varargin)

arch = computer('arch');

switch lower(arch)
    case 'maci64' %%% PERSONAL MAC:
        homepath = '/Users/neil/';
        matlabpath = '/Users/neil/Drive/MATLAB/';
    case 'glnxa64' %%% EEPC-291 or BALENA
        cd('~/');
        p = strsplit(pwd,'/');
        if any(strcmpi(p,'home')) % BALENA
            homepath = '/home/f/nh351/';
            matlabpath = '/home/f/nh351/MATLAB/';
        else % EEPC-291
            homepath = '/u/f/nh351/';
            matlabpath = '/u/f/nh351/Documents/MATLAB/';
        end
end

if ~isempty(varargin)
    if strcmpi(varargin{:},'matlab')
        path = matlabpath;
    else
        path = homepath;
    end
else
    path = homepath;
end













