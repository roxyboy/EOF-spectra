% this script takes the ensemble u, v output 
% from mitgcm runs and creates averaged fields of 
% two point correlations.  
addpath /tank/chaocean/MITgcm/utils/matlab;
clear;
file_ocean='/tank/chaocean/qjamet/RUNS/ORAR/memb';
file_out='/tank/chaocean/bill/RUNS/ORAR/ensmean/';
yearfirst=1970;yearlast=1995;
yearyearfirst=1963;
DRF=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/DRF');
DRC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/DRC');
RAC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/RAC');
nx=length(RAC(:,1));ny=length(RAC(1,:));nz=length(DRF);
memberfirst=0;
memberlast=35;
nomem=length(memberfirst:memberlast);
firstfirstit=790560;
%target data first iteration in 1967
%firstit is the iteration count for the first record in year yearfirst

%compute mean u, v
uum=zeros(nx*ny*nz*2,1);
for member=memberfirst:memberlast;
		if member < 10
			membo=['0' int2str(member)];
		else
			membo=int2str(member);
		end
	filein=[file_ocean membo '/run1967/ocn/diag_ocnTave.0001421280.data'];
fid=fopen(filein,'r','b');
fseek(fid,2*1000*900*46*4,'bof');
uu=fread(fid,[2000*900*46],'real*4');
uum=uum+uu;
member
%next member
end;

uum=uum/nomem;

fid=fopen([file_out 'utestjan1967'],'w','b');
fwrite(fid,uum,'real*4');

