
% nph_text.m

% Valid uses:
% nph_text([X Y],'string',varargin)
% nph_text([X Y W H],'string',varargin)
% nph_text(gca,[X Y],'string',varargin)
% nph_text(gcf,[X Y W H],'string',varargin)

% Units are always normalized to the size of the figure/axes on the paper.

% EDIT: OMG I ADDED A BORDER!!!!
% nph_text(...,'textborder','k',optionalextraspec)

function OUT = nph_text(varargin)

% n = 1;

% check for handle:
if strcmpi(class(varargin{1}),'matlab.graphics.axis.Axes')
    h = varargin{1};
    varargin = varargin(2:end);
else
    h = gca; % default to axes rather than figure.
end

% expect numeric position vector:
switch length(varargin{1})
    case 2
        pos = [varargin{1} 0.1 0.1];
        fitbox = 'On';
        horzalign = 'center';
        vertalign = 'middle';
    case 4
        pos = varargin{1};
        fitbox = 'Off';
        horzalign = 'left';
        vertalign = 'bottom';
    otherwise
        error('Either 2- or 4-element position vectors supported. Units are normalised to axes size and position')
end

% Expect string next:
if ischar(varargin{2})
    str = varargin{2};
else
    str = ' ';
end

% See if user added anything:
if length(varargin) > 2
    extraspec = varargin(3:end);
else
    extraspec = {};
end

% Check for textborder:
if any(strcmpi(extraspec,'textborder'))
    borderind = find(strcmpi(extraspec,'textborder'));
    borderval = extraspec{borderind+1};
    borderflag = 1;
    l = ones(size(extraspec));
    l(borderind) = 0;
    l(borderind+1) = 0;
    extraspec = extraspec(logical(l));
else
    borderflag = 0;
end

% Calculate position:
apos = h.Position;
tpos = [apos(1)+(pos(1)*apos(3)) apos(2)+(pos(2)*apos(4)) apos(3:4).*pos(3:4)];

% Text Default Specifications:
textspec = {'Units','normalized','linestyle','none','backgroundcolor',[1 1 1],'facealpha',0,'horizontalalignment',horzalign,'verticalalignment',vertalign,'FitBoxToText',fitbox};

% Try to add border... (method copied from textborder.m)
if borderflag
    %border around the text, composed of 4 text objects
    offsets = [0 -1; -1 0; 0 1; 1 0];
    Borders = struct;
    for k = 1:4
        % h = text(x, y, string, 'Color',border_color, varargin{:});
        Borders.(['h' num2str(k)]) = annotation('textbox','position',tpos,'string',str,textspec{:},extraspec{:},'color',borderval);
        Borders.(['h' num2str(k)]).Units = 'pixels';
        Borders.(['h' num2str(k)]).Position(1:2) = Borders.(['h' num2str(k)]).Position(1:2) + offsets(k,:);
        
%         %add offset in pixels
%         set(h, 'Units','pixels')
%         ppos = get(h, 'Position');
%         set(h, 'Position', [ppos(1:2) + offsets(k,:), 0 0])
%         set(h, 'Units','normalized')

    end
    
    %then the actual text inside the border
    OUT = annotation('textbox','position',tpos,'string',str,textspec{:},extraspec{:});

%     %same process as above but with 0 offset; corrects small roundoff
%     %errors
%     set(OUT, 'Units','pixels')
%     ppos = get(OUT, 'Position');
%     set(OUT, 'Position', [ppos(1:2), 0 0])
%     set(OUT, 'Units','normalized')
    
else
    % Normal Text:
    OUT = annotation('textbox','position',tpos,'string',str,textspec{:},extraspec{:});
end

















