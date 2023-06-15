close all
clear all

%define units
mm=1e-3;
um=1e-6;
nm=1e-9;
c=3e8;

% initialize pump parameters and functions
lp=.775e-6; %pump wavelength
wp = 2*pi*c/lp; 
kp=2*pi*ktpy(lp/um)/lp;
lam = @(w) 2*pi*c./w;
% convert to FWHM
pbi =(2/(2*sqrt(2*log(2))))*.6e-9; % pump bandwidth
pb =(2*pi*c*pbi)/lp^2; %convert to freq
aw =@(ws,wi) exp(-((ws+wi-wp)/(pb)).^2); %pump envelope function

%initialize crystal parameters and function
%KTP refractive index from SPDcalc
ktpx = @ (x) sqrt(2.10468 + 0.89342*x.^2./(x.^2-0.04438)-0.01036.*x.^2);
ktpy = @(x) sqrt(2.0993 + 0.922683*x.^2./(x.^2-0.0467695)-0.0138408.*x.^2);
ktpz = @(x) sqrt(1.9446 + 1.3617*x.^2./(x.^2-0.047)-0.01491.* x.^2);
k=@(w,n) 2*pi*n(lam(w)/um)./lam(w);
dk = @(ws,wi) k(ws+wi,ktpy) - k(ws,ktpy) - k(wi,ktpz); %phase matching
Lcryst=27.38*mm; % crystal length
LamC=abs(2*pi/dk(wp/2,wp/2)); %poling period

% Create wavelength grid around degenerate 1550nm
Nlam=100; % wavelength points
dLam=10; % wavelength extent
LamAll = [(1550-dLam/2):(dLam)/(Nlam-1):(1550+dLam/2)]*nm;
wall=2*pi*c./LamAll;
[Ws,Wi] = meshgrid(wall,flip(wall));

%Create crystal in spatial domain and take FFT
%FFT and sampling rate have to be/ multiples of each other
N=50; % samples per poling period
LN=20; % crystal lengths in entire spatial window (pads with zeros)
Lall=Lcryst*LN;
Fs = 1/(LamC/N);            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = (Lall)*Fs;             % Length of signal
t = (0:L-1)*T;        % Time vector
polL=LamC; % poling period

t=t-Lall/2; % shift crystal position to zero

% Calculate maximum duty cycle based on minimum poling period
% mind=0;
mind=3.425e-6;
mind=mind/4; % not sure why this is required
absDmax=1-mind/((LamC));

% create apodization function
FWHM =14.60*mm; % apodization width
apd = FWHM/(2*sqrt(2*log(2))); % convert to 1/e
cheiff = @(x) exp(-x.^2/(apd)^2); % effective nonlinearity


% construct poling amplitude in spatial domain
tN = @(x) (x*Fs+Lcryst*Fs/2); %points in sapce
Nw= Lcryst*Fs; % number of points 
% function which outputs a square wave by taking the sign of a sin wave
% with duty cycle determine by the effective nonlinearity
sW= @(x,f) sign(sin(x.*(2*pi*f))+duty(x,Lcryst,cheiff,absDmax));

% calculate function over spatial domain
S = sW(t,(1/polL));
% cut out function over crystal region
winC=and(t>-Lcryst/2,t<+Lcryst/2);

% convert crystal to a series of domain widths
ti=find(diff(S.*winC) ~= 0)+1;
widths=diff(ti)*T;
%round the widths to .025um (given by raicol)
widths=widths/(.025*um);
widths= round(widths)*.025*um;
crys=[];

%construct list of widths and convert back to the spatial domain
for i=1:size(widths,2)
    sign=1-2*mod(i,2);
    crys=horzcat(crys,sign*ones(1,round(widths(i)/T)));
    crys=horzcat(crys,0*ones(1,round(.025*um/T)));
end

%pad the crystal
PadN=size(t,2)-size(crys,2);
cryst = [zeros(1,PadN/2),crys,zeros(1,PadN/2)];

figure(2)
subplot(1,2,1)
plot(t,S.*winC,t,cryst) %plot both crystals, cheiff and duty cycle
hold on
hold on
plot(t,duty(t,Lcryst,cheiff,absDmax))
plot(t,cheiff(t))
legend('polling','duty','cheiff')
xlim([-Lcryst,Lcryst])


title('Crystal Poling')
xlabel('z (m)')
ylabel('chi(z)')


%Take FFT of crystal
% Y = fft(S.*winC);
Y=fft(cryst);
P2 = abs(Y/L); 
% conversion due to DFT
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
% define spatial frequencies
f = Fs*(0:(L/2))/L;

%plot phase matching function
figure(2)
subplot(1,2,2)
plot(2*pi*f,P1) 
hold on
plot(2*pi*f,max(P1).*abs(sinc(Lcryst*(2*pi*f-2*pi/LamC)/(2*pi)))) 
title('Phase Matching Function')
xlabel('dk (1/m)')
ylabel('h(dk)')
xlim([2*pi/LamC*(1-.01),2*pi/LamC*(1+.01)])


% calculate JSI
JSA2 = phiwF(Ws,Wi).*aw(Ws,Wi);
s2=svd(JSA2);
lambda2=s2.^2;
lambda2=lambda2./sum(lambda2);
K2=1/sum(lambda2.^2);
Leff = trapz(t,cheiff(t).*winC)/Lcryst;

% normalize JSA and plot
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




function [d]= duty(x,Lcryst,cheiff,absDmax)
% this function calculates the duty cycle to acheive a certain effective
% nonlinearity. The duty cycle is limited to some absolute duty cycle
% (-1,0,1)-> (0,100%)
cheiffx=cheiff(x);
d1=-sqrt(1-cheiffx);
d1(d1<-absDmax)=-absDmax;
d2=sqrt(1-cheiffx);
d2(d2>absDmax)=absDmax;
d=(x<0).*d1+(x>0).*d2;
end    

