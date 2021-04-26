

% Assigns the SGWEX Modelling box into the base workspace for convenience


function varargout = sgwexbox(varargin)

% % % SGWEX BOX:
% minlon = -42.5; maxlon = -31.71;
% minlat = -58.55; maxlat = -50.46;

% WARPING-CORRECTED SGWEX BOX:
% True centre location of box.
centrelat = -54.5269874182436;
centrelon = -37.1367495437688;

% Box width ~1200x900, diagonals = 1500km
% 4 --- 3
% |     |
% 1 --- 2

% NOTE: don't use the referenceEllipsoid/referenceSphere functions, I can't
% seem to get a correct answer with them and i don't know why. Had this
% trouble at Dstl too.

[boxlat(1),boxlon(1)] = reckon(centrelat,centrelon,km2deg(750),180+45);
[boxlat(2),boxlon(2)] = reckon(centrelat,centrelon,km2deg(750),180-45);
[boxlat(3),boxlon(3)] = reckon(centrelat,centrelon,km2deg(750),+45);
[boxlat(4),boxlon(4)] = reckon(centrelat,centrelon,km2deg(750),-45);

if nargout == 1
    
    sgwexbox.centrelon = centrelon;
    sgwexbox.centrelat = centrelat;
    sgwexbox.boxlon = boxlon;
    sgwexbox.boxlat = boxlat;
    
    varargout{1} = sgwexbox;
    
else
    
    assignin('base','centrelon',centrelon);
    assignin('base','centrelat',centrelat);
    
    assignin('base','boxlon',boxlon);
    assignin('base','boxlat',boxlat);
    
end


% Plot?
if nargin > 0
    if any(strcmpi('plot',varargin{:}))
        
        hold on; grid on; whitefig;
        plot(boxlon([1:end 1]),boxlat([1:end 1]),'color','r','linewi',2);
        
    end
end

