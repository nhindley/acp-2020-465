

% Saturate an RGB color triplet or Nx3 colourmap by the factor FACTOR.

% EDIT: Also adjusted to work on images

function OUT = nph_saturate(IN,factor)

IN_hsv = double(rgb2hsv(IN));

sz = size(IN);

switch sum(sz ~= 1)
    case 1
        % single color
        IN_hsv(2) = factor * IN_hsv(2);
        IN_hsv(IN_hsv > 1) = 1;
    case 2
        % Nx3 colormap
        IN_hsv(:,2) = factor * IN_hsv(:,2);
        IN_hsv(IN_hsv > 1) = 1;
    case 3
        % MxNx3 image
        IN_hsv(:,:,2) = factor * IN_hsv(:,:,2);
        
end

OUT = hsv2rgb(IN_hsv);

switch class(IN)
    case 'uint8'
        OUT = 255 * OUT;
end

% ensure class is presevered:
OUT = cast(OUT,'like',IN);









