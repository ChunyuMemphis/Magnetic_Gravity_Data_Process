

close all
clear all
up = 1; %%% upward 1 km 
% INPUT THE INCLINATION in degrees HERE----------------------
I_deg=64;
%-----------------------------------------------
% AND THE DECLINATION in degrees HERE----------------------
D_deg=1.1;
%----------------------------------------

high_limit=80000; %plotting values for corrected data
low_limit=-60000;


I=(I_deg*pi)/180;
D=(D_deg*pi)/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   Original data, but sampled %%%%%%%%%%%%%%%%%%%%%%%%%%%
A=load('magnetic_anomaly.csv');           
diur=A(1:20:end,:); 

i=sqrt(-1);

nx=128;
ny=nx;

linx(1,:)=linspace(min(diur(:,1)),max(diur(:,1)),nx);
liny(:,1)=linspace(min(diur(:,2)),max(diur(:,2)),ny);

IBA=griddata(diur(:,1),diur(:,2),diur(:,3),linx,liny,'v4');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   Detrend %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IBA1 = detrend_2d(IBA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Reduce to pole %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx=mean(diff(linx));
dy=mean(diff(liny));

nyqx=(1/(2*dx));
nyqy=(1/(2*dy));

kx=linspace(-nyqx,nyqx,nx);
ky=linspace(-nyqy,nyqy,ny);

L=cos(I)*cos(D);
l=cos(I)*cos(D);
R=cos(I)*sin(D);
r=cos(I)*sin(D);
Q=sin(I);
q=sin(I);

[KX KY]=meshgrid(kx,ky);
Kx=KX.*KX;
Ky=KY.*KY;

for m= 1:nx
    for n= 1:ny
        RTP(m,n)=(Kx(m,n)+Ky(m,n))/(( (i*L*KX(m,n)) + (i*R*KY(m,n)) + (Q*((Kx(m,n)+Ky(m,n)).^0.5))) .* ( (i*l*KX(m,n)) + (i*r*KY(m,n)) + (q*((Kx(m,n)+Ky(m,n)).^0.5)) ));
    end
end

RTP=fftshift(RTP);

filtR_dat=real(ifft2(fft2(IBA).*RTP));%%based on Gunn 1975 method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filt_up=zeros(nx,ny);
filt_down=zeros(nx,ny);
filt_derivative=zeros(nx,ny);
x_derivative=zeros(nx,ny);
y_derivative=zeros(nx,ny);
n=1;%%% the nth derivative
for i=1:nx
    for j=1:ny
        filt_up(i,j)=exp((-1*up).*(((kx(1,i).^2)+(ky(1,j).^2)).^0.5));
        filt_down(i,j)=exp(up.*(((kx(1,i).^2)+(ky(1,j).^2)).^0.5));
        filt_derivative(i,j)=kx(1,i).^2+ky(1,j).^2;
        x_derivative(i,j)=(i*kx(1,i))^n;
        y_derivative(i,j)=(i*ky(1,i))^n;
    end
end

g_up=real(ifft2(fft2(IBA).*fftshift(filt_up)));
g_down=real(ifft2(fft2(IBA).*fftshift(filt_down)));
derivative=real(ifft2(fft2(IBA).*fftshift(filt_derivative)));
x_deriv=real(ifft2(fft2(IBA).*fftshift(x_derivative)));
y_deriv=real(ifft2(fft2(IBA).*fftshift(y_derivative)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:128
     linx(i,:)=linx(1,:);
end
long=reshape(linx,[16384,1]);
liny=liny';
 for i=2:128
     liny(i,:)=liny(1,:);
 end
liny=liny';
lat=reshape(liny,[16384,1]);

M1=reshape(IBA,[16384,1]);
data1=[long,lat,M1];
dlmwrite('original.txt',data1,'delimiter','\t','precision','%.4f')

M2=reshape(IBA1,[16384,1]);
data2=[long,lat,M2];
dlmwrite('detrend.txt',data2,'delimiter','\t','precision','%.4f')

M3=reshape(filtR_dat,[16384,1]);
data3=[long,lat,M3];
dlmwrite('RTP.txt',data3,'delimiter','\t','precision','%.4f')

M4=reshape(g_up,[16384,1]);
up=[long,lat,M4];
dlmwrite('up.txt',up,'delimiter','\t','precision','%.4f')

M5=reshape(g_down,[16384,1]);
down=[long,lat,M5];
dlmwrite('down.txt',down,'delimiter','\t','precision','%.4f')

M6=reshape(derivative,[16384,1]);
vertical=[long,lat,M6];
dlmwrite('vertical.txt',vertical,'delimiter','\t','precision','%.4f')

M7=reshape(x_deriv,[16384,1]);
xderivative=[long,lat,M7];
dlmwrite('xderivative.txt',xderivative,'delimiter','\t','precision','%.4f')

M8=reshape(y_deriv,[16384,1]);
yderivative=[long,lat,M8];
dlmwrite('yderivative.txt',yderivative,'delimiter','\t','precision','%.4f')







