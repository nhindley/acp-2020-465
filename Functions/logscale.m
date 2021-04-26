

% set log scale for current axes, depending on x,y,z axis.

% logscale('x');
% logscale('x','y');
% logscale('x','y','z');


function logscale(varargin)

switch nargin
    case 0
        set(gca,'xscale','log')
    otherwise
        for i = 1:length(varargin)
            ax = varargin{i};
            switch ax
                case 'x'
                    set(gca,'xscale','log')
                case 'y'
                    set(gca,'yscale','log')
                case 'z'
                    set(gca,'zscale','log')
            end
        end
end




