

% Draws coastline into current axes using 10m coastlines shapefile.

% INPUTS:
%
%  - BoundingBox should have the format:
% [minlon minlat ; maxlon maxlat]
% which defines the box within which to consider coastlines
%
% - varargin can be anything you would put into a plot() command


% nph_draw_coastline([LonBox(1) LatBox(1) ; LonBox(2) LatBox(2)],AltOffset,MinIslandSize,coastlinespec{:});




function OUT = nph_draw_coastline(BoundingBox,AltOffset,varargin)

% DETERMINE ARCHITECTURE:
direc = localpath('matlab');

if ~isempty(varargin)
    if isnumeric(varargin{1})
        MinIslandSize = varargin{1};
        varargin = varargin(2:end);
    else
        MinIslandSize = 0;
    end
    
%     if any(strcmpi(varargin,'mplot'))
%         mplotflag = 1;
%         varargin(strcmpi(varargin,'mplot')) = [];
%     else
%         mplotflag = 0;
%     end
    
    if any(strcmpi(varargin,'noplot'))
        plotflag = 0;
        varargin(strcmpi(varargin,'noplot')) = [];
    else
        plotflag = 1;
    end
    
else
    MinIslandSize = 0;
%     mplotflag = 0;
    plotflag = 1;
end

% ShapeFile = [direc 'topography/ne/ne_10m_coastline/ne_10m_coastline.shp'];
ShapeFile = 'ne_10m_coastline.shp'; % as long as this is somewhere in the matlab path it'll find it!

coastline = shaperead(ShapeFile,'BoundingBox',double(BoundingBox),'UseGeoCoords',true);

if plotflag
    
    for i=1:length(coastline)
        
        if length(coastline(i).Lon) < MinIslandSize
            continue
        end
        
        altvec = AltOffset.*ones(size(coastline(i).Lon));
        
%         if mplotflag
%             hold on; m_plot(coastline(i).Lon,coastline(i).Lat,varargin{:},'clipping','on');
%         else
            hold on; plot3(coastline(i).Lon,coastline(i).Lat,altvec,varargin{:},'clipping','on');
%         end
        
        drawnow;
        
    end
    
else
    
    OUT = coastline;
    
end






