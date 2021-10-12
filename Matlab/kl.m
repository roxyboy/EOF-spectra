% this script takes the ensemble u, v output 
% from mitgcm runs and creates averaged fields of 
% two point correlations.  
addpath /tank/chaocean/MITgcm/utils/matlab;
clear;
file_ocean='/tank/chaocean/qjamet/RUNS/ORAR/memb';
file2_ocean='/tank/chaocean/bill/RUNS/ORAR/poje/';
DRF=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/DRF');
RAC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/RAC');
RC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/RC');
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
%xfirst=351;xlast=398;
%yfirst=668;ylast=715;
xfirst=300;xlast=xfirst+47;
yfirst=649;ylast=yfirst+47;
xr=xfirst:xlast;
yr=yfirst:ylast;
%which depth?
nr=10;
xl=length(xr);yl=length(yr);
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
uums=uc(xr,yr);
vvms=vc(xr,yr);

%two point calcs for  u, v
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
fclose(fid);
uu=reshape(uu,nx,ny,nz,2);
u=uu(:,:,:,1);
v=uu(:,:,:,2);
%move to mass points
%select depth nr 
uc(1:nx-1,:)=(u(2:nx,:,nr)+u(1:nx-1,:,nr))/2;uc(nx,:)=uc(1,:);
vc(:,1:ny-1)=(v(:,2:ny,nr)+v(:,1:ny-1,nr))/2;vc(:,ny)=vc(:,1);
%subselect region
us=uc(xr,yr);vs=vc(xr,yr);
%remove mean, weight according to area
%us=(us-uums).*sqrt(RACs);
%vs=(vs-vvms).*sqrt(RACs);
us=us-uums;
vs=vs-vvms;
	us=reshape(us,xl*yl,1);
	vs=reshape(vs,xl*yl,1);
	usl=[us vs];
	x(:,member-memberfirst+1)=reshape(usl,2*xl*yl,1);
%next member
member
end
cov_mat = x(:,1:membernumber-1)*x(:,1:membernumber-1)';
[V,D] = eigs(cov_mat,membernumber-1);         %eigen values of cov matrix

area=sum(sum(RACs,2),1);
D=D/(membernumber-1)/area/2;
figure(1);
clf;
plot(diag(D));
set(gca,'xlim',[1 membernumber-1]);
xlabel('Mode Number');
ylabel('KE (m^2/s^2)');
dth=int2str(abs(floor(RC(nr))));
title('Modal Spectrum');
saveas(gcf,['/tank/chaocean/bill/RUNS/ORAR/ensmean/modalspectrum_35dof_d_' dth '.pdf']);
savefig(['/tank/chaocean/bill/RUNS/ORAR/ensmean/modalspectrum_35dof_d_' dth '.fig']);
fid=fopen(['/tank/chaocean/bill/RUNS/ORAR/ensmean/eigfjan1967_35dof_d_' dth],'w','b');
fwrite(fid,V,'real*4');
fclose(fid);
fid=fopen(['/tank/chaocean/bill/RUNS/ORAR/ensmean/eigvjan1967_35dof_d_' dth],'w','b');
fwrite(fid,D,'real*4');
fclose(fid);
figure(1);
clf;
subplot(3,2,1);
eu=reshape(V(1:xl*yl,end),xl,yl);
ev=reshape(V(xl*yl+1:2*xl*yl,end),xl,yl);
longitude=360-98+(0:1300)/12;
latitude=-20+(0:899)/12;
cs=contour(longitude(xr),latitude(yr),eu',20);colorbar;
axis square;
set(gca,'xtick',10:2:360);
set(gca,'xticklabel',-([10:2:360]-360));
ylabel('Latitude (N)');
title(['u POD 1, D=' dth]);
subplot(3,2,2);
cs=contour(longitude(xr),latitude(yr),ev',20);colorbar;
axis square;
set(gca,'xtick',10:2:360);
set(gca,'xticklabel',-([10:2:360]-360));
title(['v POD 1, D=' dth]);
subplot(3,2,3);
eu=reshape(V(1:xl*yl,end-1),xl,yl);
ev=reshape(V(xl*yl+1:2*xl*yl,end-1),xl,yl);
cs=contour(longitude(xr),latitude(yr),eu',20);colorbar;
axis square;
set(gca,'xtick',10:2:360);
set(gca,'xticklabel',-([10:2:360]-360));
ylabel('Latitude (N)');
title(['u POD 2']);
subplot(3,2,4);
cs=contour(longitude(xr),latitude(yr),ev',20);colorbar;
axis square;
title('v POD 2');
subplot(3,2,5);
eu=reshape(V(1:xl*yl,end-2),xl,yl);
ev=reshape(V(xl*yl+1:2*xl*yl,end-2),xl,yl);
cs=contour(longitude(xr),latitude(yr),eu',20);colorbar;
axis square;
set(gca,'xtick',10:2:360);
set(gca,'xticklabel',-([10:2:360]-360));
xlabel('Longitude (W)');
ylabel('Latitude (N)');
title('u POD 3');
subplot(3,2,6);
cs=contour(longitude(xr),latitude(yr),ev',20);colorbar;
set(gca,'xtick',10:2:360);
set(gca,'xticklabel',-([10:2:360]-360));
axis square;
xlabel('Longitude (W)');
title('v POD 3');
saveas(gcf,['/tank/chaocean/bill/RUNS/ORAR/ensmean/pods1-3_35dof_d_' dth '.pdf']);
savefig(['/tank/chaocean/bill/RUNS/ORAR/ensmean/pods1-3_35dof_d_' dth '.fig']);

