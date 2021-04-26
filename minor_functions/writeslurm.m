
% Slurm = struct;
% Slurm.Filename = '~/Desktop/writeslurmtest.slurm';
% 
% Slurm.Name = 'fart';
% Slurm.Batch = 'batch-short';
% 
% Slurm.Command = {'module load matlab/2016b','airs_3dst_z40km_BALENA',matlab -nodisplay -nosplash -r "PRO_airs_3dst_bias_binning; exit"};

% writes a slurm file to the filename in the Slurm structure:

function writeslurm(Slurm)

f = fieldnames(Slurm);

if ~any(strcmpi(f,'Filename'))
    Slurm.Fileame = [pwd '/Slurm.slurm'];
end
if ~any(strcmpi(f,'Name'))
    Slurm.Name = 'Name';
end
if ~any(strcmpi(f,'Batch'))
    Slurm.Batch = 'batch-all';
end
if ~any(strcmpi(f,'Command'))
    Slurm.Command = 'disp(''No command found'')';
end
if ~any(strcmpi(f,'Batch'))
    Slurm.Batch = 'batch-all';
end
if ~any(strcmpi(f,'WallTime'))
    Slurm.WallTime = '360'; % wall time in minutes
end

% now write slurm:

fid = fopen(Slurm.Filename,'w');

writeline(fid,'#!/bin/bash')
nph_newline(fid);

writeline(fid,'## Name of the job')
writeline(fid,['#SBATCH --job-name=' Slurm.Name])
nph_newline(fid);

writeline(fid,'## Account to charge to')
writeline(fid,'#SBATCH --account=free')
nph_newline(fid);

writeline(fid,'## Walltime required')
writeline(fid,['#SBATCH --time=' num2str(Slurm.WallTime) ':00'])
nph_newline(fid);

writeline(fid,'## Number of node required and tasks per node')
writeline(fid,'#SBATCH --nodes=1')
writeline(fid,'#SBATCH --ntasks-per-node=16')
nph_newline(fid);

writeline(fid,'## Where to direct the output to')
writeline(fid,['#SBATCH --output=log.' Slurm.Name '.%j.out'])
writeline(fid,['#SBATCH --error=log.' Slurm.Name '.%j.err'])
writeline(fid,['#SBATCH --partition=' Slurm.Batch])
nph_newline(fid);

% Process commands:
writeline(fid,'## Load Matlab environment')
writeline(fid,'module load matlab/2016b')
writeline(fid,['matlab -nodisplay -nosplash -r "' Slurm.Command '; exit"'])
nph_newline(fid);

fclose(fid);

end






















