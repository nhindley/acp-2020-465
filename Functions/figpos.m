
% set figure position on screen

% IN    - if 4-element vector [x y width height], place on screen and interpret as pixels.
%       - if 2-element vector [width height], centre in middle of screen and make height
%       width the inputs as normalized (inputs between 0 and 1)

% EDIT 2020: Now scaling everything to the aspect ratio of a 1920x1200
% monitor (yeah, 16:10 baby!)


function figpos(IN)

asp = 1920./1200; % 1.6

switch isnumeric(IN) % allow for just screensize
    case 0
        if any(strcmpi(IN,{'screensize','fullscreen'}))
            set(gcf,'position',get(0,'screensize'));
        end
    case 1 % or numeric wavlues specifying widths and heights
        switch length(IN)
            case 2
                scrn = get(0,'screensize');
                current_asp = scrn(3)./scrn(4);
                modifier = asp./current_asp;
                xwid = IN(1) * scrn(3) * modifier;
                ywid = IN(2) * scrn(4);
                x = (scrn(3)/2) - (xwid/2);
                y = (scrn(4)/2) - (ywid/2);
                set(gcf,'position',[x y xwid ywid]);
            case 4
                set(gcf,'position',IN);
        end
        
end
end