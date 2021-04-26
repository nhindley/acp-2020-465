
% can't remember whether to use fliplr or flipud for 1D column/row vectors?
% Use reverse!

function OUT = reverse(IN)

% % get orientation:
sz = size(IN);

switch length(sz(sz ~= 1))
    case 1 % 1D
        OUT = IN(end:-1:1);
    case 2 % 2D
        IN = squeeze(IN);
        OUT = IN(end:-1:1,end:-1:1);
    case 3 % 3D
        IN = squeeze(IN);
        OUT = IN(end:-1:1,end:-1:1,end:-1:1);
end

% And finally for scalar input, do nothing:
if length(IN) == 1
    OUT = IN;
end

% put back into original orientation:
% OUT = reshape(OUT,sz);








