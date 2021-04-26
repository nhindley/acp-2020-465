

% writes a new line in a text file
% not to be confused with newline.m, a built in matlab one which just
% generates the character sprintf('\n').

function nph_newline(fileID,varargin)

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


% % % 
% % % % newline:
% % % arch = computer;
% % % switch upper(arch)
% % %     case {'PCWIN','PCWIN64'}
% % %         fprintf(fileID,'\r\n');
% % %     case {'GLNX86','GLNXA64','MACI64'}
% % %         fprintf(fileID,'\n');
% % %     otherwise
% % %         fprintf(fileID,'\n');
% % % end

end








