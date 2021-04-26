   

% for a given altitude, what's the vertical resolution of Lars's 3D AIRS
% retrieval?

% these are estimated average values!


function OUT = airs3d_vert_res(varargin)

switch nargin
    case 0
        z = [0:3:60 65 70];
        interpmethod = 'linear';
    case 1
        if isnumeric(varargin{1})
            z = varargin{1};
            interpmethod = 'linear';
        else
            z = [0:3:60 65 70];
            interpmethod = varargin{1};
        end
    case 2
        z = varargin{1};
        interpmethod = varargin{2};
end

alt = [0:3:60 65 70];

res = [
        0  ,25
        3  ,25
        6  ,25
        9  ,25
        12 ,20
        15 ,15
        18 ,10.25;...
        21 ,6.5;...
        24 ,6.75;...
        27 ,7.25;...
        30 ,7.5;...
        33 ,7.75;...
        36 ,8.75;...
        39 ,9;...
        42 ,9;...
        45 ,10;...
        48 ,10;...
        51 ,11;...
        54 ,13.5;...
        57 ,14.5;...
        60 ,13.25;...
        65 ,13.5;...
        70 ,19.5];
    
F = griddedInterpolant(alt,res(:,2),interpmethod,'nearest');

OUT = F(z);





















