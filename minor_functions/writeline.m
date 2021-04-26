

% Writes a string to text file and makes a new line


function writeline(fileID,string,varargin)

fprintf(fileID,'%s',string);

% NEW LINE!
if ~isempty(varargin)
    arch = varargin{1};
else
    arch = 'mac';
end

switch lower(arch)
    case {'win','windows'}
        fprintf(fileID,'\r\n');
    case {'mac','apple'}
        fprintf(fileID,'\n');
    case {'linux','unix','gnu'}
        fprintf(fileID,'\n');
    otherwise
        fprintf(fileID,'\n');
end


end






























