
% Use a Tukey-style tapered window around a chosen range:
% (specified in elements atm...)

% We expect IN to be a 3D matrix of XYZ dimensions, and we apply the window
% to the THIRD dimension.

% You can have indeces that extend off the top of the matrix, but not of
% the bottom, i.e. indices greater than length() but not < 1.

function OUT = sgwex_airs_apply_height_window(IN,start_index,end_index,tapering_length)

if nargin == 3
    tapering_length = 6;
end

if isodd(tapering_length)
    tapering_length = tapering_length + 1;
end

sz = size(IN);

% number of elements equal to one:
n_eq2one = (end_index - start_index) - tapering_length;

% create window section
tkwin = nph_tukeywin(tapering_length,n_eq2one,tapering_length);

n_leading_zeros = (start_index-tapering_length/2);

n_trailing_zeros = size(IN,3) - (length(tkwin) + n_leading_zeros);

winvec = [zeros(n_leading_zeros,1) ; tkwin ; zeros(n_trailing_zeros,1)] ;

% apply
win = permute(repmat(winvec,1,sz(1),sz(2)),[2 3 1]);

% trim to region of input.
win = win(1:sz(1),1:sz(2),1:sz(3));

OUT = IN .* win; 






















