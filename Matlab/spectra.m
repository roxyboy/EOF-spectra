clear
figure
clf;
addpath /tank/chaocean/MITgcm/utils/matlab;
for t=1:2;
t
if t == 1
depth = 94;
else
depth=628;
end;
ext=['_35dof_d_' int2str(depth)];
pl1=-5/3;
pl2=-3;
fid=fopen(['/tank/chaocean/bill/RUNS/ORAR/ensmean/lengthsjan1967' ext],'r','b');
L=fread(fid,'real*4');
fclose(fid);
fid=fopen(['/tank/chaocean/bill/RUNS/ORAR/ensmean/waveksjan1967' ext],'r','b');
kll=fread(fid,'real*4');
kx=kll(1:35);ly=kll(36:end);
fclose(fid);
%sort wavenumbers and convert to cyc/km
[ks,ikx]=sort(kx(1:end-1));
[lsort,jky]=sort(ly(1:end-1));
ks=ks*1.e3;lsort=lsort*1.e3;
fid=fopen(['/tank/chaocean/bill/RUNS/ORAR/ensmean/eigvjan1967' ext],'r','b');
ev=fread(fid,'real*4');
fclose(fid);
ev=reshape(ev,35,35);
evv=diag(ev);
evv=evv(1:end-1);
subplot(2,2,2*(t-1)+1);
plot(log10(ks(1:end)),log10(evv(ikx(1:end))),'-x');
hold on;
plot([log10(ks(3)) log10(ks(30))],[log10(evv(3)) log10(evv(3))+pl1*(log10(ks(30))-log10(ks(3)))],'r');
plot([log10(ks(3)) log10(ks(30))],[log10(evv(3)) log10(evv(3))+pl2*(log10(ks(30))-log10(ks(3)))],'k');
set(gca,'xlim',[-3 -2]);
set(gca,'xtick',[-3:.2:-2]);
set(gca,'xticklabel',[1000 630 400 250 160 100]);
xlabel('cycles/km');
ylabel('Energy (J)');
title(['Zonal Spectra, ' int2str(depth) ' m']);
text(-2.9,-11,'-3','Color','k','fontsize',14);
text(-2.9,-11.5,'-5/3','Color','r','fontsize',14);
grid on;
hold off;
subplot(2,2,2*(t-1)+2);
plot(log10(lsort(1:end)),log10(evv(jky(1:end))),'-x');
hold on;
plot([log10(lsort(3)) log10(lsort(30))],[log10(evv(3)) log10(evv(3))+pl1*(log10(lsort(30))-log10(lsort(3)))],'r');
plot([log10(lsort(3)) log10(lsort(30))],[log10(evv(3)) log10(evv(3))+pl2*(log10(lsort(30))-log10(lsort(3)))],'k');
text(-2.8,-11,'-3','Color','k','fontsize',14);
text(-2.8,-11.5,'-5/3','Color','r','fontsize',14);
set(gca,'xlim',[-3 -2]);
set(gca,'xtick',[-3:.2:-2]);
set(gca,'xticklabel',[1000 630 400 250 160 100]);
xlabel('cycles/km');
ylabel('Energy (J)');
title(['Meridional Spectra, ' int2str(depth) ' m']);
grid on;
hold off;
t
end;
savefig(['/tank/chaocean/bill/RUNS/ORAR/ensmean/allspectra']);
saveas(gca,['/tank/chaocean/bill/RUNS/ORAR/ensmean/allspectra'],'pdf');

