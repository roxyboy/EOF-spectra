clear
addpath ../scripts
DXC=rdmds('../data/DXC');
DYC=rdmds('../data/DYC');
XC=rdmds('../data/XC');
YC=rdmds('../data/YC');

%read mean u, v
file=['../data/ens/diag_ocnTave1421280'];
fid=fopen(file,'r','b');
fseek(fid,2*1000*900*46*4,'bof');
uu=fread(fid,[2000*900*46],'real*4');fclose(fid);
uu=reshape(uu,1000,900,46,2);
u=uu(:,:,:,1);v=uu(:,:,:,2);clear uu;
%define region
xfirst=300;xlast=xfirst+47;
yfirst=649;ylast=yfirst+47;
xr=xfirst:xlast;
yr=yfirst:ylast;
xl=length(xr);
yl=length(yr);

%cosine factor
cosf=cos(2*3.14159*YC(xr(10),yr)/360);

%depths
nr=[10 21];
%compute eddy fields
memberfirst=0;memberlast=35;
membernumber=length(memberfirst:memberlast);
filename='../data/diag_ocnTave';
mh=0;
for member=memberfirst:memberlast;
    mh=mh+1;mh
    if member < 10
        membo=['0' int2str(member)];
    else
        membo=int2str(member);
    end;
    file=[filename membo '.0001421280.data'];
    fid=fopen(file,'r','b');
    fseek(fid,2*1000*900*46*4,'bof');
uu=fread(fid,[2000*900*46],'real*4');fclose(fid);
uu=reshape(uu,1000,900,46,2);
um=uu(:,:,:,1);vm=uu(:,:,:,2);clear uu;
um=um-u;vm=vm-v;
for kz=1:2;
ums(:,:,kz,mh)=um(xr,yr,nr(kz));vms(:,:,kz,mh)=vm(xr,yr,nr(kz));
%weight by cosf;
%for j=1:length(yr);
%    umsc(:,j,kz)=ums(:,j,kz,mh)*sqrt(cosf(j));
%    vmsc(:,j,kz)=vms(:,j,kz,mh)*sqrt(cosf(j));
%end;
%detrend zonally
%slope
umsx=(ums(end,:,kz,mh)-ums(1,:,kz,mh))/(xl-1);
vmsx=(vms(end,:,kz,mh)-vms(1,:,kz,mh))/(xl-1);
for j=1:yl;
    umsdx(:,j)=ums(:,j,kz,mh)-umsx(j)*(xr-xr(1))';
     vmsdx(:,j)=vms(:,j,kz,mh)-vmsx(j)*(xr-xr(1))';
end;
%detrend meridionally
umsy=(umsdx(:,end)-umsdx(:,1))/(yl-1);
vmsy=(vmsdx(:,end)-vmsdx(:,1))/(yl-1);
for i=1:xl;
    umsdxdy(i,:,kz,mh)=umsdx(i,:)-umsy(i)*(yr-yr(1));
    vmsdxdy(i,:,kz,mh)=vmsdx(i,:)-vmsy(i)*(yr-yr(1));
end;
    coru(:,:,kz,mh)=ums(:,:,kz,mh)-umsdxdy(:,:,kz,mh);
    corv(:,:,kz,mh)=vms(:,:,kz,mh)-vmsdxdy(:,:,kz,mh);
end;

end;   
fid=fopen('../data/udet','w','b');
fwrite(fid,umsdxdy,'real*4');fclose(fid);
fid=fopen('../data/vdet','w','b');
fwrite(fid,vmsdxdy,'real*4');fclose(fid);
fid=fopen('../data/ucor','w','b');
fwrite(fid,coru,'real*4');fclose(fid);
fid=fopen('../data/vcor','w','b');
fwrite(fid,corv,'real*4');fclose(fid);
fid=fopen('../data/ueddy','w','b');
fwrite(fid,ums,'real*4');fclose(fid);
fid=fopen('../data/veddy','w','b');
fwrite(fid,vms,'real*4');fclose(fid);

for kz=1:2;
for mh=1:membernumber;
umsh=fft2(ums(:,:,kz,mh));
umshs=fftshift(umsh);
vmsh=fft2(vms(:,:,kz,mh));
vmshs=fftshift(vmsh);
totsp(:,:,kz,mh)=umshs.*conj(umshs)+vmshs.*conj(vmshs);
umsdxdyh=fft2(umsdxdy(:,:,kz,mh));
umsdxdyhs=fftshift(umsdxdyh);
vmsdxdyh=fft2(vmsdxdy(:,:,kz,mh));
vmsdxdyhs=fftshift(vmsdxdyh);
dxdysp(:,:,kz,mh)=umsdxdyhs.*conj(umsdxdyhs)+vmsdxdyhs.*conj(vmsdxdyhs);
coruh=fft2(coru(:,:,kz,mh));
coruhs=fftshift(coruh);
corvh=fft2(corv(:,:,kz,mh));
corvhs=fftshift(corvh);
corsp(:,:,kz,mh)=corvhs.*conj(corvhs)+coruhs.*conj(coruhs);
crosssp1=coruhs.*conj(umsdxdyhs)+conj(coruhs).*umsdxdyhs;
crosssp2=corvhs.*conj(vmsdxdyhs)+conj(corvhs).*vmsdxdyhs;
crosssp(:,:,kz,mh)=crosssp1+crosssp2;
end;mh,end;

corspm=sum(corsp,4)/36;
crossspm=sum(crosssp,4)/36;
dxdyspm=sum(dxdysp,4)/36;
totspm=sum(totsp,4)/36;

