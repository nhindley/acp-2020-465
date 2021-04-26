

% 1-D cross-track polynomial fit.

% with no dimension argument specified: we expect IN to be AIRS Temperature in XTxATxAlt dimensions, that is, we will be taking
% the fit across the FIRST dimension.

% [perturbations,background] = airs_remove_bg(IN,dim,polyorder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT JULY 2020
% I've noticed that the polynomial fit seems to be sensitive to the
% direction of the data that's put in in - that is, you get a different fit
% if you put the same data in but the order of the elements is reversed.
% This shouldn't happen! Bad!
% FIXED - I think. I changed the dummy spatial vector to be balanced around
% zero.
% I also took out the mean removal step once I'd fixed this. It could be
% thrown off by weird data so didn't want to risk it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function varargout = nph_airs_4dp_detrend(varargin)


switch nargin
    case 1
        IN = varargin{1};
        dim = 1;
        polyorder = 4;
    case 2
        IN = varargin{1};
        dim = varargin{2};
        polyorder = 4;
    otherwise
        IN = varargin{1};
        dim = varargin{2};
        polyorder = varargin{3};
end

% permute to a big line of data in the chosen dimension:
sz = size(IN);
dims = 1:1:numel(sz);
order = unique([dim,dims],'stable');
IN2 = permute(IN,order);
IN2 = reshape(IN2,[sz(order(1)),prod(sz(order(2:end)))]);

BG = nan(size(IN2));

% new nph version of code below:
ssz = size(IN2);

% % % % EDIT: Now removed vvvvv
% % % % one key difference is that the mean is calculated seperately, rather than
% % % % fitted by the matrix division.
% % % % One reason this seems to be helping is that by making the mean zero, you
% % % % don't get some of the higher order terms trying to fit the mean and they
% % % % are left to fit the higher order features. This shouldn't work, but it
% % % % seems to, so perhaps I don't understand how the maths works :)


% FIRST create dummy spatial vector
% x = (1:ssz(1))';
x = linspace(-1,1,ssz(1))';


% for each cross-track row, or equivalent...
for i = 1:ssz(2)
    
    % data to be fitted
    t = IN2(:,i);
    
    % DEAL WITH NANS
    % simply exclude them from the fit data and the coordinate space
    nanlocs = isnan(t);
    if any(nanlocs)
        t_in = t(~nanlocs);
        x_in = x(~nanlocs);
    else
        t_in = t;
        x_in = x;
    end
    
    %     % find and remove the mean:
    %     tm = mean(t);
    %     t = t - tm;
    
    % make vectors of powers * spatial vector
    M = bsxfun(@power,x_in,0:polyorder);
    
    % use matrix division to estimate fit coefficients
    warning off
    coeffs = mldivide(M,t_in);
    
    % use polyval to evaluate these coefficients on the dummy grid, NaNs
    % and all:
    bg = polyval(reverse(coeffs),x);
    
    % replace NaNs in the output?
    bg(nanlocs) = NaN;
    
    %     % put the mean back on
    %     bg = bg + tm;
    
    % and assign:
    BG(:,i) = bg;
    
end


%put back into original shape
BG = reshape(BG,[sz(order(1)),(sz(order(2:end)))]);
BG = permute(BG,order);

% Remove perts to leave background:
perts = IN - BG;


% Mutliple outputs?
switch nargout
    case 1
        varargout{1} = perts;
    case 2
        varargout{1} = perts;
        varargout{2} = BG;
end
















% % % % % % % Defaults:
% % % % % % if ~isempty(varargin)
% % % % % %     if isnumeric(varargin{1})
% % % % % %         polyorder = varargin{1};
% % % % % %     end
% % % % % % else
% % % %     polyorder = 4;
% % % % % % end
% % % %
% % % % sz = size(IN);
% % % %
% % % % % first, reshape for fun:
% % % % IN = reshape(IN,[sz(1) sz(2)*sz(3)]);
% % % %
% % % % % Set up fittype for X and Y directions:
% % % % ft = fittype(['poly' num2str(polyorder)]);
% % % %
% % % % % Initialise background
% % % % perts = zeros(size(IN));
% % % %
% % % % % coord frame:
% % % % x = (1:sz(1))';
% % % %
% % % % % Find polynomial fit:
% % % % for i = 1:(sz(2)*sz(3))
% % % %
% % % %     y = IN(:,i);
% % % %
% % % % % Fit model to data.
% % % % [fitresult, ~] = fit(x,y,ft,'Normalize','on');
% % % %
% % % % % Get perts from original minus background:
% % % % perts(:,i) = y - feval(fitresult,x);
% % % %
% % % % % Give zero mean: (should be close to zero anyway)
% % % % perts(:,i) = perts(:,i) - mean(perts(:,i));
% % % %
% % % % end
