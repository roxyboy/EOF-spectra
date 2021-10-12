%%%%%%%%%%
% Plotting of the 3D EOF spectra.
%%%%%%%%%%
clear
addpath /tank/chaocean/MITgcm/utils/matlab;
fid=fopen('figs/uveofs','r','b');
uveof=fread(fid,'real*4');fclose(fid);
uveof=reshape(uveof,48,48,2,35,2);
uv=squeeze(uveof(:,:,1,:,:));
vv=squeeze(uveof(:,:,2,:,:));
fid=fopen('figs/uveigs_36','r','b');
uveigs=fread(fid,'real*4');fclose(fid);
uveigs=reshape(uveigs,35,2);
DXC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/DXC');
RC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/RC');
DYC=rdmds('/tank/chaocean/grid_chaO/gridMIT_update1/DYC');
%define x and y ranges;
xfirst=300;xlast=xfirst+47;
yfirst=649;ylast=yfirst+47;
xr=xfirst:xlast;
yr=yfirst:ylast;
dxcs=DXC(xr,yr);
dycs=DYC(xr,yr);
xl=length(xr);yl=length(yr);
%which depth?
nr=[10 21];
depth = [RC(nr(1)) RC(nr(2))];
for kz=1:length(nr);
for k=1:35;
uvx(:,:,kz)=squeeze((uv(2:end,:,k,kz)-uv(1:end-1,:,k,kz))./dxcs(1:end-1,:));
uvy(:,:,kz)=squeeze((uv(:,2:end,k,kz)-uv(:,1:end-1,k,kz))./dycs(:,1:end-1));
vvx(:,:,kz)=squeeze((vv(2:end,:,k,kz)-vv(1:end-1,:,k,kz))./dxcs(1:end-1,:));
vvy(:,:,kz)=squeeze((vv(:,2:end,k,kz)-vv(:,1:end-1,k,kz))./dycs(:,1:end-1));
k2(k,kz)=sum(sum(uvx(:,:,kz).^2+vvx(:,:,kz).^2,2),1);
l2(k,kz)=sum(sum(uvy(:,:,kz).^2+vvy(:,:,kz).^2,2),1);
end;
end;
Lx=sqrt(1./k2)*2*pi;
Ly=sqrt(1./l2)*2*pi;
L=[Lx Ly];
keff=sqrt(k2)/2/pi;
leff=sqrt(l2)/2/pi;
kkk=[keff leff];
%kkk in cycles/m
%create 3-d plot
dkeff(1:34,:)=(keff(2:35,:)-keff(1:34,:))*1.e3;
dleff(1:34,:)=(leff(2:35,:)-leff(1:34,:))*1.e3;
uveigsn=abs(uveigs(1:34,:));

for m=1:2;
	subplot(2,1,m);
	hold on;
	lgk(:,m)=log10(keff(:,m)*1e3);
	lgl(:,m)=log10(leff(:,m)*1e3);
	if m == 1
	plot3(lgk(1:34,m),lgl(1:34,m),-4*ones(34,1),'.b');
else
	plot3(lgk(1:34,m),lgl(1:34,m),-5*ones(34,1),'.b');
end;
view(110,15);
hold on;
if m == 1
	set(gca,'zlim',[-4 -.5])
else
	set(gca,'zlim',[-5 -1.5]);
end;
	set(gca,'xlim',[-3 -2]);
	set(gca,'ylim',[-2.8 -2]);
	plot3(lgk(1:33,m),lgl(1:33,m),log10(uveigsn(1:33,m)),'mx');
plot3(-3*ones(33,1),lgl(1:33,m),log10(uveigsn(1:33,m)),'r.');
plot3(lgk(1:33,m),-2.8*ones(33,1),log10(uveigsn((1:33),m)),'g.');
if m == 1
plot3([-2.6 -2.2],[-2.8 -2.8],[-1 -2.2],'k');
plot3([-3 -3],[-2.6 -2.2],[-1 -2.2],'k');
else
plot3([-2.6 -2.2],[-2.8 -2.8],[-2 -3.2],'k');
plot3([-3 -3],[-2.6 -2.2],[-2 -3.2],'k');
end;
set(gca,'xtick',[-3:.5:-2],'xticklabel',[-3:.5:-2]);
set(gca,'ytick',[-2.8:.4:-2],'yticklabel',[-2.8:.4:-2]);
xlabel('log(k)' );
ylabel('log(l)');
zlabel('log(E)');
axis square;
if m == 1
	title('94m');
else
	title('628m');
end;
grid on;
hold off;
end;
savefig('figs/eofspectra');
saveas(gcf,'figs/eofspectra.pdf');
