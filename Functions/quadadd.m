

% Add in quadrature.
% Should be array safe :)

function OUT = quadadd(varargin)

IN = varargin(:);

% initialise:
s = zeros(size(IN{1}));

% sum of squared:
for i = 1:nargin
    s = s + IN{i}.^2;
end

% take root:
OUT = sqrt(s);











