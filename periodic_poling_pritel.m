close all
clear all

mm=1e-3;
um=1e-6;

ktpx = @ (x) sqrt(2.10468 + 0.89342*x.^2./(x.^2-0.04438)-0.01036.*x.^2);
ktpy = @(x) sqrt(2.0993 + 0.922683*x.^2./(x.^2-0.0467695)-0.0138408.*x.^2);
ktpz = @(x) sqrt(1.9446 + 1.3617*x.^2./(x.^2-0.047)-0.01491.* x.^2);

c=3e8;
lp=.775e-6;
wp = 2*pi*c/lp;
kp=2*pi*ktpy(lp/um)/lp;
lam = @(w) 2*pi*c./w;
pbi =(2/(2*sqrt(2*log(2))))*.6e-9;
pb =(2*pi*c*pbi)/lp^2;
k=@(w,n) 2*pi*n(lam(w)/um)./lam(w);
aw =@(ws,wi) exp(-((ws+wi-wp)/(pb)).^2)
% aw =@(ws,wi) exp(-((ws+wi-wp)/(2*pb)).^2)
dk = @(ws,wi) k(ws+wi,ktpy) - k(ws,ktpy) - k(wi,ktpz)
% Lcryst=22.58*mm;
Lcryst=27.38*mm;
% Lcryst=.52*mm;
LamC=abs(2*pi/dk(wp/2,wp/2));
% LamC=46.52*um
phiw= @(ws,wi) sinc(Lcryst*(dk(ws,wi)+2*pi/LamC)/(2*pi)) %matlab adds a pi

% FWHM=21*mm;
% apd = FWHM/(2*sqrt(2*log(2)));
% apod = @(x) exp(-(x-Lcryst/2).^2/(sqrt(2)*apd)^2);


LamAll = [1545:(10)/(99):1555]*1e-9;
wall=2*pi*c./LamAll;
[Ws,Wi] = meshgrid(wall,flip(wall));
figure(1)
subplot(2,3,1)
imagesc(LamAll*1e6,LamAll*1e6,abs(aw(Ws,Wi)))
title('Pump Envelope')
subplot(2,3,2)
imagesc(LamAll*1e6,LamAll*1e6,abs(phiw(Ws,Wi)))
title('Phase Matching')
subplot(2,3,3)
imagesc(LamAll*1e6,LamAll*1e6,(abs(phiw(Ws,Wi)).*abs(aw(Ws,Wi))).^2)
title('JSI')
JSA1 = abs(phiw(Ws,Wi)).*aw(Ws,Wi);
s1=svd(JSA1);
lambda1=s1.^2;
lambda1=lambda1./sum(lambda1);
K1=1/sum(lambda1.^2)



%FFT and sampling rate have to be/ multiples of each other
mm=1e-3;
um=1e-6;

N=50;
Lall=Lcryst*20;
Fs = 1/(LamC/N);            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = (Lall)*Fs;             % Length of signal
t = (0:L-1)*T;        % Time vector
polL=LamC;

t=t-Lall/2;

% mind=5e-6;
% mind=10e-6;
% mind=0;
mind=3.425e-6;
mind=mind/4;
absDmax=1-mind/((LamC));

%Guassian Window
% cheiff =@(x) exp(-x.^2/(Lcryst/5)^2);
FWHM =14.60*mm;
% FWHM=Lcryst*.6;
% FWHM=Lcryst;
apd = FWHM/(2*sqrt(2*log(2)));
cheiff = @(x) exp(-x.^2/(apd)^2);



tN = @(x) (x*Fs+Lcryst*Fs/2);
Nw= Lcryst*Fs;

% cheiff = @(x) ChebChi(x,Lcryst*1.6,1);
% cheiff = @(x) tukChi(x,Lcryst*1.6,.9)
% cheiff = @(x) (1-abs((tN(x)-Nw/2)/(ratio*Nw/2)));

% sigma=.8;
% p=8;
% cheiff = @(x) exp(-((tN(x)-Nw/2)/(sigma*Nw/2)).^p);
% k=1:Nw;
% W0=@(k) 
% chebstuff =@(x) W0(k).*(-exp(1i*pi/(N+1))).^(k).*exp(1i.*2.*pi.*k.*tN(x)./(Nw+1))
% 1/(Nw+1)*sum(chebstuff)



% cheiff=@(x) exp(-0*x.^2/(Lcryst/5)^2)
% cheiff =@(x) double(and(x<Lcryst/4,x>-Lcryst/4));
% duty = @(x) asin(cheiff(x))/pi;
sW= @(x,f) sign(sin(x.*(2*pi*f))+duty(x,Lcryst,cheiff,absDmax));
% sW = @(x,f)  sign(sin(x.*(2*pi*f))).*cheiff(x);
% sW = @(x,f)  sign(sin(x.*(2*pi*f))).*.5;

S = sW(t,(1/polL));
% S = (square(2*pi*t/polL,50));
winC=and(t>-Lcryst/2,t<+Lcryst/2);


ti=find(diff(S.*winC) ~= 0)+1;
widths=diff(ti)*T;
widths=widths/(.025*um);
widths= round(widths)*.025*um;
crys=[];
for i=1:size(widths,2)
    sign=1-2*mod(i,2);
    crys=horzcat(crys,sign*ones(1,round(widths(i)/T)));
    crys=horzcat(crys,0*ones(1,round(.025*um/T)));
end
PadN=size(t,2)-size(crys,2);

cryst = [zeros(1,PadN/2),crys,zeros(1,PadN/2)];

figure(2)
subplot(1,2,1)
plot(t,S.*winC,t,cryst)
% plot(t,cryst)
hold on
hold on
plot(t,duty(t,Lcryst,cheiff,absDmax))
plot(t,cheiff(t))
legend('polling','duty','cheiff')
xlim([-Lcryst,Lcryst])


title('Crystal Poling')
xlabel('z (m)')
ylabel('chi(z)')



% Y = fft(S.*winC);
Y=fft(cryst);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure(2)
subplot(1,2,2)
plot(2*pi*f,P1) 
hold on
plot(2*pi*f,max(P1).*abs(sinc(Lcryst*(2*pi*f-2*pi/LamC)/(2*pi)))) 
title('Phase Matching Function')
xlabel('dk (1/m)')
ylabel('h(dk)')
xlim([2*pi/LamC*(1-.01),2*pi/LamC*(1+.01)])

figure(1)
subplot(2,3,4)
imagesc(LamAll*1e6,LamAll*1e6,aw(Ws,Wi))
title('Pump Envelope (Apod)')
subplot(2,3,5)
phiwF =@(ws,wi) interp1(2*pi*f,P1,-dk(ws,wi));
imagesc(LamAll*1e6,LamAll*1e6,phiwF(Ws,Wi))
title('Phase Matching (Apod)')
subplot(2,3,6)
imagesc(LamAll*1e6,LamAll*1e6,(phiwF(Ws,Wi).*aw(Ws,Wi)))
title('JSA (Apod)')

JSA2 = phiwF(Ws,Wi).*aw(Ws,Wi);
s2=svd(JSA2);
lambda2=s2.^2;
lambda2=lambda2./sum(lambda2);
K2=1/sum(lambda2.^2)
Leff = trapz(t,cheiff(t).*winC)/Lcryst

lambGs=2*pi*3e8./(Ws)-1.55*um;
lambGi=2*pi*3e8./(Wi)-1.55*um;
JSI = JSA2^2/(max(max(JSA2^2)));
JSI=JSA2/(max(max(JSA2)));

JSI=JSA2/(max(max(JSA2)));
f=figure
f.Position = [100 100 250 200]
imagesc(LamAll*1e6,LamAll*1e6,1-JSA2)
colormap bone
xlabel('Signal (\mum)')
ylabel('Idler (\mum)')
title('JSA (Apod)')


JSI=JSA1/(max(max(JSA1)));
f=figure
f.Position = [100 100 250 200]
imagesc(LamAll*1e6,LamAll*1e6,1-JSI)
colormap bone
xlabel('Signal (\mum)')
ylabel('Idler (\mum)')
title('JSA')
% surf(lambGs,lambGi,abs(JSI-.5))
% locs=abs(JSI-.5)<.01;
% fwhms=lambGs(locs);
% fwhmi=lambGi(locs);
% fwhm = sqrt(fwhms.^2+fwhmi.^2)*1e9;
% fwhmM=mean(fwhm)*2;
% transition_indices = find(diff(S.*winC) ~= 0) + 1;
% max_widths = diff(transition_indices(1:2:end)) * T;
% min_widths = diff(transition_indices(2:2:end)) * T;


% txt=fileread('export.json');
% data=jsondecode(txt);
% x=data.plots.xvals;
% y=data.plots.yvals;
% JSIw=data.plots.zvals;
% s2w=svd(sqrt(JSIw));
% lambdaw=s2w.^2;
% lambdaw=lambdaw./sum(lambdaw);
% Kw=1/sum(lambdaw.^2)

% figure
% subplot(1,3,1)
% contourf(log(JSA1.^2),10)
% subplot(1,3,2)
% contourf(log(JSA2.^2),10)
% subplot(1,3,3)
% JSIw=flip(JSIw).';
% contourf(log(JSIw),10)
% imagesc(abs(JSA2.^2/max(max(JSA2.^2))-JSIw/max(max(JSIw))))
function [d]= duty(x,Lcryst,cheiff,absDmax)
% cheiff =@(x) exp(-x.^2/(Lcryst/5)^2);
cheiffx=cheiff(x);
d1=-sqrt(1-cheiffx);
d1(d1<-absDmax)=-absDmax;
d2=sqrt(1-cheiffx);
d2(d2>absDmax)=absDmax;
d=(x<0).*d1+(x>0).*d2;

% for i = 1:max(size(x))
%     if x(i)<0 
% %     if 1==1
% %         d(i)=(absDmax/.5)*asin(cheiffx(i))/pi-absDmax;
%         d(i)=max([-sqrt(1-cheiffx(i)^2),-absDmax]);
%     else
% %         d(i)=absDmax-(absDmax/.5)*asin(cheiffx(i))/pi;
%         d(i)=min([sqrt(1-cheiffx(i)^2),absDmax]);
% %         d(i)=1-cheiffx(i);
%     end
% end
end

function [chi]=ChebChi(x,Width,dB)
dx=mean(diff(x));
Nwin=Width/dx;
win = chebwin(Nwin,100*dB);
Nwin2=max(size(win));
% xwin=[-Width/2:(Width)/(Nwin-1):Width/2];
xwin = [1:Nwin2]*dx-Width/2;
chi=interp1(xwin,win,x);
chi(isnan(chi))=0;
end

% function [cryst]=gen_from_widths(t,widths,T)
% crys=[];
% for i=1:size(widths,2)
%     sign=1-2*mod(i,2);
%     horzcat(crys,sign*ones(1,round(widths(i)/T)));
% end
% plot(crys)
% figure
% PadN=size(t,2)-size(crys,2);
% cryst = [zeros(1,PadN/2),crys,zeros(1,PadN/2)];
% end
    

function [chi]=tukChi(x,Width,dB)
dx=mean(diff(x));
Nwin=Width/dx;
win = tukeywin(Nwin,dB);
Nwin2=max(size(win));
% xwin=[-Width/2:(Width)/(Nwin-1):Width/2];
xwin = [1:Nwin2]*dx-Width/2;
chi=interp1(xwin,win,x);
chi(isnan(chi))=0;
end

