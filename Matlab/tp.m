% this script takes the ensemble u, v output 
% from mitgcm runs and creates averaged fields of 
% two point correlations.  
addpath /tank/chaocean/MITgcm/utils/matlab;
clear;
file_ocean='/tank/chaocean/qjamet/RUNS/ORAR/memb';
file2_ocean='/tank/chaocean/bill/RUNS/ORAR/poje/';
DRF=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/DRF');
RAC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/RAC');
nx=length(RAC(:,1,1));ny=length(RAC(1,:,1));
nz=length(DRF);
yearfirst=1970;yearlast=1995;
yearyearfirst=1963;
memberfirst=0;
memberlast=35;
membernumber=length(memberfirst:memberlast);
firstfirstit=790560;yearyearfirst=1963;
%target data first iteration in 1967
yearfirst=1967;
firstit=firstfirstit+(yearfirst-yearyearfirst)*73*2160;
%firstit is the iteration count for the first record in year yearfirst
%define x and y ranges;
xfirst=351;xlast=398;
yfirst=668;ylast=715;
xr=xfirst:xlast;
yr=yfirst:ylast;
%which depth?
nr=10;
xl=length(xr);yl=length(yr);
tpc=zeros(xl*yl*2,xl*yl*2);
%get mean velocities
uu=zeros(nx*ny*nz*2,1);
fid=fopen('/tank/chaocean/bill/RUNS/ORAR/ensmean/utestjan1967','r','b');
uu=fread(fid,'real*4');
uu=reshape(uu,nx,ny,nz,2);

%move to mass points
%select depth nr 
uc(1:nx-1,:)=(uu(2:nx,:,nr,1)+uu(1:nx-1,:,nr,1))/2;uc(nx,:)=uc(1,:);
vc(:,1:ny-1)=(uu(:,2:ny,nr,2)+uu(:,1:ny-1,nr,2))/2;vc(:,ny)=vc(:,1);

%subselect data
	RACs=RAC(xr,yr);
	RACs=reshape(RACs,xl*yl,1);
	RACs=[RACs RACs];
	RACs=reshape(RACs,2*xl*yl,1);
uums=uc(xr,yr);
vvms=vc(xr,yr);

%two point correlations for  u, v
for member=memberfirst:memberlast;
		if member < 10
			membo=['0' int2str(member)];
		else
			membo=int2str(member);
		end
filein=[file_ocean membo '/run1967/ocn/diag_ocnTave.000' int2str(firstit) '.data'];
fid=fopen(filein,'r','b');
fseek(fid,2*1000*900*46*4,'bof');
uu=fread(fid,[2000*900*46],'real*4');
uu=reshape(uu,nx,ny,nz,2);
u=uu(:,:,:,1);
v=uu(:,:,:,2);
%move to mass points
%select depth nr 
uc(1:nx-1,:)=(u(2:nx,:,nr)+u(1:nx-1,:,nr))/2;uc(nx,:)=uc(1,:);
vc(:,1:ny-1)=(v(:,2:ny,nr)+v(:,1:ny-1,nr))/2;vc(:,ny)=vc(:,1);
us=uc(xr,yr);vs=vc(xr,yr);
us=us-uums;
vs=vs-vvms;
datau(:,:,member-memberfirst+1)=us;
datav(:,:,member-memberfirst+1)=vs;

%next member
member
end
datau=datau.*RACs;
datav=datav.*RACs;
datau_map=eof(datau);
datav_map=eof(datav);
