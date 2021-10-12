%%%%%%%%%%
% This script takes the ensemble u, v output
% from mitgcm runs and computes the eigenvalues/eigenvectors of the
% two point correlations.
% Requires the MITgcm Matlab toolkit to function.
%%%%%%%%%%
addpath /tank/chaocean/MITgcm/utils/matlab;
clear;
file_ocean='/tank/chaocean/qjamet/RUNS/ORAR/memb';
file2_ocean='/tank/chaocean/bill/RUNS/ORAR/ensmean/';
DRF=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/DRF');
RAC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/RAC');
RC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/RC');
XC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/XC');
YC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/YC');
DYC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/DYC');
DXC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/DXC');
DYG=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/DYG');
DXG=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/DXG');
nx=length(RAC(:,1,1));ny=length(RAC(1,:,1));
nz=length(DRF);
memberfirst=0;
memberlast=35;
memvec=[memberfirst:memberlast];
membernumber=length(memvec);
firstfirstit=790560;yearyearfirst=1963;
%target data first iteration in 1967
yearfirst=1967;
firstit=firstfirstit+(yearfirst-yearyearfirst)*73*2160;
%firstit is the iteration count for the first record in year yearfirst
xfirst=300;xlast=xfirst+47;
yfirst=649;ylast=yfirst+47;
xr=xfirst:xlast;
yr=yfirst:ylast;
rac=RAC(xr,yr);
cosf=cos(2*3.14159*YC(xr(10),yr)/360);
%compute domain area
area=sum(sum(rac,2),1);
% depths
nr=[10 21];
xl=length(xr);yl=length(yr);
%get mean velocities
uu=zeros(nx*ny*nz*2,1);
fid=fopen(['/tank/chaocean/bill/RUNS/ORAR/ensmean/run' int2str(yearfirst) '_' int2str(membernumber) '/oce/diag_ocnTave' int2str(firstit)],'r','b');
fseek(fid,2*1000*900*46*4,'bof');
uu=fread(fid,[2000*900*46],'real*4');fclose(fid);
uu=reshape(uu,nx,ny,nz,2);

%move to mass points
%select depth nr
uc(1:nx-1,:,:)=(uu(2:nx,:,:,1)+uu(1:nx-1,:,:,1))/2;uc(nx,:,:)=uc(1,:,:);
vc(:,1:ny-1,:)=(uu(:,2:ny,:,2)+uu(:,1:ny-1,:,2))/2;vc(:,ny,:)=vc(:,1,:);

%subselect data
uums=uc(xr,yr,:);
vvms=vc(xr,yr,:);

%earth radium in m
a=6370e3;

%angular resolution in lon, lat
dlon=1/12*3.14159/180;
dlat=1/11.62*3.14159/180;

%two point calcs for  u, v
kenergy=zeros(length(nr),1);
kenergy1=kenergy;
%for member=memberfirst:memberlast;
imem=0;
for member=memvec;
imem=imem+1;
		if member < 10
			membo=['0' int2str(member)];
		else
			membo=int2str(member);
		end
filein=[file_ocean membo '/run1967/ocn/diag_ocnTave.000' int2str(firstit) '.data'];
fid=fopen(filein,'r','b');
fseek(fid,2*1000*900*46*4,'bof');
uu=fread(fid,[2000*900*46],'real*4');
fclose(fid);
uu=reshape(uu,nx,ny,nz,2);
u=uu(:,:,:,1);
v=uu(:,:,:,2);
%move to mass points
uc(1:nx-1,:,:)=(u(2:nx,:,:)+u(1:nx-1,:,:))/2;uc(nx,:,:)=uc(1,:,:);
vc(:,1:ny-1,:)=(v(:,2:ny,:)+v(:,1:ny-1,:))/2;vc(:,ny,:)=vc(:,1,:);
%subselect region
us=uc(xr,yr,:);vs=vc(xr,yr,:);
us=us-uums;
vs=vs-vvms;
%compute ke
for kz=1:2;
kenergy1(kz)=kenergy1(kz)+sum(sum((us(:,:,nr(kz)).^2+vs(:,:,nr(kz)).^2).*DXC(xr,yr).*DYC(xr,yr),2),1)/2;
end;
	%weight by sqrt cos
	for kz=1:length(nr);
	for i=1:length(xr);
		us(i,:,nr(kz))=us(i,:,nr(kz)).*sqrt(cosf);
		vs(i,:,nr(kz))=vs(i,:,nr(kz)).*sqrt(cosf);
	end;
end;
%regional energy
for kz=1:length(nr);
kenergy(kz)=kenergy(kz)+sum(sum(us(:,:,nr(kz)).^2 + vs(:,:,nr(kz)).^2,2),1)*a*a*dlon*dlat/2;
end;
	for kz=1:length(nr);
	usa=reshape(us(:,:,nr(kz)),xl*yl,1);
	vsa=reshape(vs(:,:,nr(kz)),xl*yl,1);
	usl=[usa vsa];
	x(:,imem,kz)=reshape(usl,2*xl*yl,1);
end;
%next member
member
end
kenergy=kenergy/(membernumber);
kenergy1=kenergy1/(membernumber);
for kz=1:length(nr);
cov_mat = x(:,1:membernumber-1,kz)*x(:,1:membernumber-1,kz)'*a*a*dlat*dlon/(membernumber-1)/2;
[V,D] = eigs(cov_mat,membernumber-1);         %eigen values of cov matrix
vh(:,:,kz)=V;dh(:,:,kz)=D;
end;
% values in dh have units of area integrated energy
% correct for domain area to end up in units of energy
dh=dh/area;
% eigenmodes in vh are orthonormalized when summed
vh=reshape(vh,length(xr),length(yr),2,membernumber-1,2);
for kz=1:membernumber-1;
dhs(kz,:)=dh(kz,kz,:);
end;
xl=XC(xr(1):xr(end),1);yl=YC(xr(1),yr(1):yr(end));
xl=360-xl;
for mode=1:2;
uhsurf=vh(:,:,1,mode,1);
vhsurf=vh(:,:,2,mode,1);
uhdeep=vh(:,:,1,mode,2);
vhdeep=vh(:,:,2,mode,2);
subplot(2,2,1);
cs=contour(xl,yl,uhsurf',20);colorbar;
set(gca,'xtick',68:72,'xticklabel',68:72);
set(gca,'ytick',35:38,'yticklabel',35:38);
title(['U Mode ' int2str(mode) ', 94m']);
ylabel('lat N');
subplot(2,2,2);
cs=contour(xl,yl,vhsurf',20);colorbar;
set(gca,'xtick',68:72,'xticklabel',68:72);
set(gca,'ytick',35:38,'yticklabel',35:38);
title(['V Mode ' int2str(mode) ', 94m']);
subplot(2,2,3);
cs=contour(xl,yl,uhdeep',20);colorbar;
set(gca,'xtick',68:72,'xticklabel',68:72);
set(gca,'ytick',35:38,'yticklabel',35:38);
title(['U Mode ' int2str(mode) ', 628m']);
ylabel('lat N');
xlabel('lon W');
subplot(2,2,4);
cs=contour(xl,yl,vhdeep',20);colorbar;
set(gca,'xtick',68:72,'xticklabel',68:72);
set(gca,'ytick',35:38,'yticklabel',35:38);
title(['V Mode ' int2str(mode) ', 628m']);
xlabel('lon W');
file=['figs/uveof' int2str(mode)];
savefig([file '.fig']);
saveas(gcf,[file '.pdf']);
end;
fid=fopen(['figs/uveofs_' int2str(membernumber)],'w','b');
fwrite(fid,vh,'real*4');fclose(fid);
fid=fopen(['figs/uveigs_' int2str(membernumber)],'w','b');
fwrite(fid,dhs,'real*4');fclose(fid);
