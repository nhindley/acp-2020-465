

% Replacement for nph_export_fig.m

% Same inputs as the built-in function saveas.m

function nph_saveas(varargin)

% No scaling applied now:
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% deal with different axes colors
set(gcf,'InvertHardCopy','Off');

% scalefactor = 5;
% pos = get(gcf,'Position');
% 
% figx = px2in(pos(3));
% figy = px2in(pos(4));
% 
% set(gcf,'PaperSize',scalefactor*[figx figy]);
% set(gcf,'PaperPosition',scalefactor*[0 0 figx figy]);

% applystyle;

saveas(varargin{:})

end

function inches = px2in(pixels)
    % 1920px ~ 20.5 inches. Measured with a high-tech method (ruler).
    inches = (pixels / 1920) * 20.5;

end




% % % Apply a style using hgexport:

function applystyle

    style = struct();
    style.Version = '1';
    style.Format = 'eps';
    style.Preview = 'none';
    style.Width = 'auto';
    style.Height = 'auto';
    style.Units = 'centimeters';
    style.Color = 'rgb';
    style.Background = 'w';          % '' = no change; 'w' = white background
    style.FixedFontSize = '20';
    style.ScaledFontSize = 'auto';
    style.FontMode = 'auto';
    style.FontSizeMin = '8'; % commented out as seemed to rescale upon execution
    style.FixedLineWidth = '1';
    style.ScaledLineWidth = 'auto';
    style.LineMode = 'auto';
    style.LineWidthMin = '0.5';
    style.FontName = 'auto';
    style.FontWeight = 'auto';
    style.FontAngle = 'auto';
    style.FontEncoding = 'latin1';
    style.PSLevel = '2';
    style.Renderer = 'auto';
    style.Resolution = 'auto';
    style.LineStyleMap = 'none';
    style.ApplyStyle = '0';
    style.Bounds = 'loose';
    style.LockAxes = 'on';
    style.ShowUI = 'on';
    style.SeparateText = 'off';

    hgexport(gcf,'temp_dummy',style,'applystyle', true);

end
%








