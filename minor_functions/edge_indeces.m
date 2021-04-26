

% Find the linear indeces of the edge rows and columns of the 2D input
% matrix. Counting clockwise.

% Saves a lot of faff!

function varargout = edge_indeces(IN,varargin)

sz = size(IN);

switch length(sz)
    case 2 % 2D input
        
        i1 = [1:sz(1) repmat(sz(1),1,sz(2)) sz(1):-1:1 ones(1,sz(2))];
        i2 = [ones(1,sz(1)) 1:sz(2) repmat(sz(2),1,sz(1)) sz(2):-1:1];
        
        if nargin > 1
            if any(strcmpi(varargin{:},'plot'))
                figure; hold all;
                IN(sub2ind(sz,i1,i2)) = 0;
                pcolor(IN);
                shading flat;
            end
        end
        
        inds = sub2ind(sz,i1,i2);
        
    case 3 % 3D input
        
        inds = [];
        
        i1 = [1:sz(1) repmat(sz(1),1,sz(2)) sz(1):-1:1 ones(1,sz(2))];
        i2 = [ones(1,sz(1)) 1:sz(2) repmat(sz(2),1,sz(1)) sz(2):-1:1];
        
        for z = 1:sz(3)
            
            inds = [inds sub2ind(sz,i1,i2,z*ones(size(i1)))];
            
        end
        
end

switch nargout
    case 1
        varargout{1} = inds;
    case 2
        varargout{1} = i1;
        varargout{2} = i2;
end



















