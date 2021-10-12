clear
file2_ocean='/tank/chaocean/bill/RUNS/ORAR/ensmean/';
fid=fopen([file2_ocean 'eigfjan1967_35dof_d_94'],'r','b');
fcts=fread(fid,'real*4');

fcts=reshape(fcts,48*48*2,35);
uv=reshape(fcts(1:48*48,:),48,48,35);
vv=reshape(fcts(48*48+1:end,:),48,48,35);

for k=1:35;
uva=uv(:,:,k)-sum(sum(uv(:,:,k),2),1)/48/48;
ftuh=fft2(uva);
[ii,jj]=find(abs(ftuh)==max(max(abs(ftuh))));
iii(:,k)=ii;jjj(:,k)=jj;
ftuhh=zeros(48,48);
ftuhh(ii,jj)=ftuh(ii,jj);
ftu(:,:,k)=ifft2(ftuhh);
[r(:,:,k),p(:,:,k)]=corrcoef(ftu(:,:,k),uv(:,:,k));
ftut(:,k)=reshape(ftu(:,:,k),48*48,1);
ftut(:,k)=ftut(:,k)-mean(ftut(:,k));
ftut(:,k)=ftut(:,k)/std(ftut(:,k));
uvt(:,k)=reshape(uv(:,:,k),48*48,1);
uvt(:,k)=uvt(:,k)-mean(uvt(:,k));
uvt(:,k)=uvt(:,k)/std(uvt(:,k));
uvr(:,:,k)=uv(:,:,k)-ftu(:,:,k);
end;
plot(squeeze(r(1,2,:).^2),'-x');
set(gca,'ylim',[0 1]);
xlabel('Mode number');
ylabel('Variance');
savefig('/tank/chaocean/bill/RUNS/ORAR/ensmean/rsquared');
saveas(gca,'/tank/chaocean/bill/RUNS/ORAR/ensmean/rsquared','pdf');
