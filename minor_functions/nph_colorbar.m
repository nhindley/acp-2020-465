

% Ever noticed how when you plot a colorbar, the axes automatically resize
% themselves, and then if you move the colorbar, they resize themselves
% back to their original position, and the colorbar sits infuriatingly over
% the axes? Ever noticed that? Well?
%
% Well I have, and it bloody annoying.
%
% So annoying in fact that i have written nph_colorbar, which does exactly
% what colorbar does but puts the colorbar neatly at the same relative
% position and size for each axis.
%
% Only works for location eastoutside I'm afraid.
%

function cbar = nph_colorbar(varargin)


% get axes positon FIRST
apos = get(gca,'Position');

% initialise colorbar
cbar = colorbar(varargin{:});
drawnow;

% poke it to reset the axes:
cbar.Position = cbar.Position;

% now set the colorbar in a nice place:
% x and y
cbar.Position(1) = apos(1) + apos(3) + 0.025*apos(3);
cbar.Position(2) = apos(2) + 0.2*apos(4);

% width and height
cbar.Position(3) = 0.04*apos(3);
cbar.Position(4) = 0.6*apos(4);



end


















