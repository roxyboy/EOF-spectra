%%%%%%%%%%
% Code to derive the effective wavenumbers
%%%%%%%%%%
clear
addpath /tank/chaocean/MITgcm/utils/matlab;
depth = 94;
ext=['_35dof_d_' int2str(depth)];
fid=fopen(['/tank/chaocean/bill/scripts/poje/figs/uveofs'],'r','b');
DXC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/DXC');
DYC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/DYC');
%define x and y ranges;
xfirst=300;xlast=xfirst+47;
yfirst=649;ylast=yfirst+47;
xr=xfirst:xlast;
yr=yfirst:ylast;
%which depth?
dxcs=DXC(xr,yr);
dycs=DYC(xr,yr);
nr=21;
xl=length(xr);yl=length(yr);
V=fread(fid,'real*4');
fclose(fid);
V=reshape(V,xl*yl*2,35,2);
uv=reshape(V(1:xl*yl,:,:),xl,yl,35,2);
vv=reshape(V(xl*yl+1:end,:,:),xl,yl,35,2);
for kz=1:2;
for k=1:35;
uvx=(uv(2:end,:,k,kz)-uv(1:end-1,:,k,kz))./dxcs(1:end-1,:);
uvy=(uv(:,2:end,k,kz)-uv(:,1:end-1,k,kz))./dycs(:,1:end-1);
vvx=(vv(2:end,:,k,kz)-vv(1:end-1,:,k,kz))./dxcs(1:end-1,:);
vvy=(vv(:,2:end,k,kz)-vv(:,1:end-1,k,kz))./dycs(:,1:end-1);;
k2(k,kz)=sum(sum(uvx.*uvx+vvx.*vvx,2),1);
l2(k,kz)=sum(sum(uvy.*uvy+vvy.*vvy,2),1);
end;
end;
Lx=sqrt(1./k2)*2*pi;
Ly=sqrt(1./l2)*2*pi;
L=[Lx Ly];
keff=sqrt(k2)/2/pi;
leff=sqrt(l2)/2/pi;
kkk=[keff leff];
%kkk in cycles/m
fid=fopen(['/tank/chaocean/bill/RUNS/ORAR/ensmean/lengthsjan1967' ext],'w','b');
fwrite(fid,L,'real*4');
fclose(fid);
fid=fopen(['/tank/chaocean/bill/RUNS/ORAR/ensmean/waveksjan1967' ext],'w','b');
fwrite(fid,kkk,'real*4');
fclose(fid);
