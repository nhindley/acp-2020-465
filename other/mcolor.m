
%
% Ever wondered what the RGBs of Matlab's default colors are for line
% plotting? Well now you can find out with mcolor.m!
%
% For a demo of the default colors, type mcolor with no arguments.
%
% For an RGB triplet, use mcolor('blue') or mcolor(1) etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OUT = mcolor(IN)

colornames = {...
    'blue'
    'orange'
    'yellow'
    'purple'
    'green'
    'light_blue'
    'red'};

colors = [...
    0   114 189 ; % blue
    217 83  25  ; % orange
    237 177 32  ; % yellow
    126 47  142 ; % purple
    119 172 48  ; % green
    77  190 238 ; % light_blue
    162 20  47 ]; % red
    
colors = colors ./ 255;

% EXAMPLE DEMO
if nargin == 0 || any(strcmpi(IN,{'demo','plot','example'}))
    figure; set(gcf,'color','w')
    for i = 1:7
        hold on; l = plot(1:10,(1:10)+i,'color',colors(i,:),'linewi',3);
        disp([sprintf('%4d',round(l.Color * 255)) '  ' colornames{i}])
    end
    return
end

% fix if IN is greater than the number of colors, just wrap to size(colors,1)...
if isnumeric(IN)
    positiveInput = (IN > 0);
    IN = mod(IN, size(colors,1));
    IN((IN == 0) & positiveInput) = size(colors,1);
    % ^ stolen from wrapTo360
end

switch isnumeric(IN)
    case 1
        if IN == 0
            OUT = [0 0 0];
        else
            OUT = colors(IN,:);
        end
    case 0
        OUT = colors(strcmpi(colornames,IN),:);
end











