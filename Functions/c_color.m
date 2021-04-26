

% outputs an RGB colour numbers for nice colours :)



function c = c_color(color)

% Reds
if strcmp(color,'red'),    c = [204 0 0]  ./255; end
if strcmp(color,'dark_red'),    c = 0.5.* [204 0 0]  ./255; end

% Oranges
if strcmp(color,'orange'), c = [217 83 25]./255; end
if strcmp(color,'dark_orange'), c = 0.5.*[255 128 0]./255; end

% Yellows
if strcmp(color,'yellow'),  c = [255 255 0] ./255; end
if strcmp(color,'mustard'),  c = [255 217 47] ./255; end
if strcmp(color,'gold'),  c = [254 215 0] ./255; end

% Greens
if strcmp(color,'light_green'),  c = [167 215 91] ./255; end
if strcmp(color,'green'),  c = [51 152 0] ./255; end
if strcmp(color,'dark_green'),  c = [51 152 0] ./255; end

% Blues
if strcmp(color,'very_light_blue'),  c = [200 210 255] ./255; end
if strcmp(color,'light_blue'),  c = [135 206 250] ./255; end
if strcmp(color,'blue'),   c = [0 114 189]./255; end
if strcmp(color,'dark_blue'),  c = 0.5.* [0 102 204] ./255; end

% Pinks and Purples
if strcmp(color,'pink'),  c = [255 182 193] ./255; end
if strcmp(color,'purple'), c = [76 0 153] ./255; end
if any(strcmp(color,{'dark_purple','deep_purple'})), c = 0.5* [76 0 153] ./255; end

% Browns
if strcmp(color,'brown'),  c = [84 48 5] ./255; end

% Whites
if strcmp(color,'white'),  c = [255 255 255] ./255; end

% Blacks and Greys
if strcmp(color,'very_light_grey'),   c = 2 * [96 96 96] ./255; end
if strcmp(color,'light_grey'),   c = 1.5 * [96 96 96] ./255; end
if strcmp(color,'grey'),   c = [96 96 96] ./255; end
if strcmp(color,'dark_grey'),   c = 0.5 * [96 96 96] ./255; end
if strcmp(color,'black'),  c = [0 0 0] ./255; end












if strcmp(color,'white'),  c = [255 255 255] ./255; end



end