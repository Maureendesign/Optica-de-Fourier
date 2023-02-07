%%
% 
% <<FILENAME.PNG>>
% 
%Resolución 2 um/pixel
close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 0.- Propiedades Generales %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = .550; %um
f = 3962000; %um
N=512;

te = 20;       %Tiempo de espocisión seg
D = 2;         %
rho = 100000; %um
sigma = 10;    %radianes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 1.- Definición de la apertura AP %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xap = (-N/2:1:N/2)*(2.128);
yap = (-N/2:1:N/2)*(2.128);
AP = zeros(length(xap));

for i=1:length(xap)
for j=1:length(yap)
if sqrt(xap(i)^2 + yap(j)^2) <= 305
AP(i,j) = 1;
end
end
end

for i=1:length(xap)
for j=1:length(yap)
if sqrt(xap(i)^2 + yap(j)^2) <= 209
AP(i,j) = 0;
end
end
end

for i=1:length(xap)
if abs(xap(i)) <= 10
AP(i,:) = 0;
end
end

for j=1:length(yap)
if abs(yap(j)) <= 10
AP(:,j) = 0;
end
end

figure; set(gca,'FontSize',14);
imagesc(AP);
pbaspect([1 1 1])
title('Apertura del telescopio','FontSize',14);
xlabel('Eje X (mm)')
ylabel('Eje Y (mm)')
saveas(gcf,'1-1_Apertura','png');

FAP = abs(fftshift(fft2(AP)));

%VECTORES DEL ESPACIO MTF
u = (-N/2:1:N/2)*(1/((N+1)*(2)));
v = (-N/2:1:N/2)*(1/((N+1)*(2)));

%VECTORES DEL ESPACIO PSF
xpsf = (-N/2:1:N/2)*2;
ypsf = (-N/2:1:N/2)*2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% 2.- Movimiento Lateral %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vimg = 0.1; %Velocidad = 0.1 um/seg
AreaRecorrida = vimg*te;
VelEf = zeros(length(xpsf));

for i=1:length(xpsf)
if abs(xpsf(i)) <= AreaRecorrida
VelEf(:,i) = xpsf(i)/AreaRecorrida;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 3.- Variaciones de Montura %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Var=zeros(length(xpsf));

for i=1:length(xpsf)
for j=1:length(ypsf)
Var(i,j) = real(1/(pi*sqrt((D*1.01)^2 - (xpsf(i)^2 + ypsf(j)^2))));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 4.- Distorsiones de la atmosfera %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(u)
for j=1:length(v)
MTF4(i,j) = exp(-(sigma^2)*(1 - exp(-((lambda*f/rho)^2)*(u(i)^2 + v(j)^2))));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% 5.- Detector pixeleado %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 1024;

PSF5 = zeros(length(xpsf));

for i=1:length(xpsf)
for j=1:length(ypsf)
PSF5(j,i) = xpsf(i) + ypsf(j);
end
end
PSF5 = rot90(PSF5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MODIFICADORES TOTALES %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSF1 = FAP;
PSF2 = VelEf;
PSF3 = Var;
MTF4 = MTF4;
PSF5 = PSF5;

MTF1 = abs(fftshift(fft2(PSF1)));
MTF2 = abs(fftshift(fft2(PSF2)));
MTF3 = abs(fftshift(fft2(PSF3)));
PSF4 = abs(fftshift(ifft2(MTF4)));
MTF5 = abs(fftshift(fft2(PSF5)));

MTF1 = MTF1/max(max(MTF1));
MTF2 = MTF2/max(max(MTF2));
MTF3 = MTF3/max(max(MTF3));
MTF4 = MTF4/max(max(MTF4));
MTF5 = MTF5/max(max(MTF5));

PSF1 = (PSF1/max(max(PSF1)));
PSF2 = (PSF2/max(max(PSF2)));
PSF3 = (PSF3/max(max(PSF3)));
PSF4 = (PSF4/max(max(PSF4)));
PSF5 = (PSF5/max(max(PSF5)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% GRÁFICAS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; set(gca,'FontSize',14);
imagesc(PSF1(200:312,200:312));
pbaspect([1 1 1])
title('PSF1 - Apertura del telescopio','FontSize',14);
xlabel('x-values (um)')
ylabel('y-values (um)')
saveas(gcf,'1-2_PSF','png');

figure; set(gca,'FontSize',14);
plot(xpsf(200:312),PSF1(200:312,257));
pbaspect([1 1 1])
title('PSF1 - Apertura del telescopio Perfil','FontSize',14);
xlabel('x-values (um)')
ylabel('PSF 1')
saveas(gcf,'1-3_PSF Perfil','png');

figure; set(gca,'FontSize',14);
surf(xpsf,ypsf,PSF2)
pbaspect([1 1 1])
%colorbar
title('PSF2 - Movimiento Lateral','FontSize',14);
xlabel('x-values')
ylabel('y-values')
zlabel('PSF 2')
saveas(gcf,'2-1_PSF','png');

figure; set(gca,'FontSize',14);
plot(xpsf,PSF2(257,:))
pbaspect([1 1 1])
%colorbar
title('PSF2 - Movimiento Lateral Perfil','FontSize',14);
xlabel('x-values')
ylabel('PSF 2')
saveas(gcf,'2-2_PSF Perfil','png');

figure; set(gca,'FontSize',14);
imagesc(PSF3)
pbaspect([1 1 1])
title('PSF3 - Variaciones de Montura','FontSize',14);
xlabel('x-values')
ylabel('y-values')
saveas(gcf,'3-1_PSF','png');

figure; set(gca,'FontSize',14);
plot(xpsf,PSF3((length(xpsf)-1)/2,:))
pbaspect([1 1 1])
%colorbar
title('PSF3 - Variaciones de Montura Perfil','FontSize',14);
xlabel('x-values')
ylabel('PSF 3')
saveas(gcf,'3-2_PSF Perfil','png');

figure; set(gca,'FontSize',14);
imagesc(MTF4)
pbaspect([1 1 1])
%colorbar
title('MTF4 - Distorsiones de la atmósfera','FontSize',14);
xlabel('u-values')
ylabel('v-values')
saveas(gcf,'4-1_MTF','png');

figure; set(gca,'FontSize',14);
plot(u,MTF4(:,257))
pbaspect([1 1 1])
%colorbar
title('MTF4 - Distorsiones de la atmósfera perfil','FontSize',14);
xlabel('u-values')
ylabel('MTF 4')
saveas(gcf,'4-2_MTF Perfil','png');

figure; set(gca,'FontSize',14);
surf(xpsf,ypsf,PSF5)
pbaspect([1 1 1])
%colorbar
title('PSF5 - Efecto del detector','FontSize',14);
xlabel('x-values')
ylabel('y-values')
zlabel('PSF 5')
saveas(gcf,'5-1_PSF','png');

figure; set(gca,'FontSize',14);
plot(xpsf,PSF5)
pbaspect([1 1 1])
%colorbar
title('PSF5 - Perfil','FontSize',14);
xlabel('x-values')
ylabel('PSF 5')
saveas(gcf,'5-2_PSF Perfil','png');

figure; set(gca,'FontSize',14);
imagesc(MTF1);
pbaspect([1 1 1])
title('MTF1 - Apertura del telescopio','FontSize',14);
xlabel('x-values (um)')
ylabel('y-values (um)')
saveas(gcf,'1-4_MTF','png');

figure; set(gca,'FontSize',14);
plot(u,abs(MTF1(257,:)));
pbaspect([1 1 1])
title('MTF1 - Apertura del telescopio Perfil','FontSize',14);
xlabel('u-values (um)')
ylabel('MTF 1')
saveas(gcf,'1-5_MTF Perfil','png');

figure; set(gca,'FontSize',14);
surf(u,v,MTF2)
pbaspect([1 1 1])
%colorbar
title('MTF2 - Movimiento Lateral','FontSize',14);
xlabel('u-values')
ylabel('v-values')
zlabel('MTF 2')
saveas(gcf,'2-3_MTF','png');

figure; set(gca,'FontSize',14);
plot(u,MTF2(257,:))
pbaspect([1 1 1])
%colorbar
title('MTF2 - Movimiento Lateral Perfil','FontSize',14);
xlabel('u-values')
ylabel('MTF 2')
saveas(gcf,'2-4_MTF Perfil','png');

figure; set(gca,'FontSize',14);
imagesc(abs(MTF3))
pbaspect([1 1 1])
title('MTF3 - Variaciones de Montura','FontSize',14);
xlabel('u-values')
ylabel('v-values')
saveas(gcf,'3-3_MTF','png');

figure; set(gca,'FontSize',14);
plot(u,MTF3(257,:))
pbaspect([1 1 1])
%colorbar
title('MTF3 - Variaciones de Montura Perfil','FontSize',14);
xlabel('u-values')
ylabel('MTF 3')
saveas(gcf,'3-4_MTF Perfil','png');

figure; set(gca,'FontSize',14);
imagesc(PSF4(200:312,200:312))
pbaspect([1 1 1])
%colorbar
title('PSF4 - Distorsiones de la atmósfera','FontSize',14);
xlabel('x-values')
ylabel('y-values')
saveas(gcf,'4-3_PSF','png');

figure; set(gca,'FontSize',14);
plot(xpsf(200:312),PSF4(257,200:312))
pbaspect([1 1 1])
%colorbar
title('PSF4 - Distorsiones de la atmósfera perfil','FontSize',14);
xlabel('x-values')
ylabel('PSF 4')
saveas(gcf,'4-4_PSF Perfil','png');

figure; set(gca,'FontSize',14);
surf(u,v,MTF5)
pbaspect([1 1 1])
%colorbar
title('MTF5 - Efecto del detector','FontSize',14);
xlabel('u-values')
ylabel('v-values')
zlabel('MTF 5')
saveas(gcf,'5-3_MTF','png');

figure; set(gca,'FontSize',14);
plot(u(200:312),abs(MTF5(257,200:312)))
pbaspect([1 1 1])
%colorbar
title('MTF5 - Perfil','FontSize',14);
xlabel('u-values')
ylabel('MTF 5')
saveas(gcf,'5-4_MTF Perfil','png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% PRODUCTOS FINALES %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1 = conv2(PSF2,PSF5);
A1 = imresize(A1,0.5);
%A1 = A1(257:769,257:769);
A1 = A1/max(max(A1));
B1 = conv2(PSF1,A1);
B1 = imresize(B1,0.5);
%B1 = B1(257:769,257:769);
B1 = B1/max(max(B1));
B2 = conv2(PSF3,PSF4);
B2 = imresize(B2,0.5);
%B2 = B2(257:769,257:769);
B2 = B2/max(max(B2));
PSFTotal = conv2(B1,B2);
PSFTotal = imresize(PSFTotal,0.5);
%PSFTotal = PSFTotal(257:769,257:769);
PSFTotal = PSFTotal/max(max(PSFTotal));

figure; set(gca,'FontSize',14);
plot(xpsf,abs(PSFTotal(:,257)))
pbaspect([1 1 1])
%colorbar
title('PSF Total en Y','FontSize',14);
xlabel('y-values')
ylabel('PSF Total')
saveas(gcf,'6-1_PSF Total Y','png');

figure; set(gca,'FontSize',14);
plot(xpsf,abs(PSFTotal(257,:)))
pbaspect([1 1 1])
%colorbar
title('PSF Total en X','FontSize',14);
xlabel('x-values')
ylabel('PSF Total')
saveas(gcf,'6-2_PSF Total X','png');

figure; set(gca,'FontSize',14);
imagesc(abs(PSFTotal))
pbaspect([1 1 1])
colorbar
title('PSF Total','FontSize',14);
xlabel('x-values')
ylabel('y-values')
saveas(gcf,'6-3_PSF Total','png');

MTFTotal = MTF1.*MTF2.*MTF3.*MTF4.*MTF5;
MTFTotal = MTFTotal/max(max(MTFTotal));

figure; set(gca,'FontSize',14);
imagesc(abs(MTFTotal(200:312,200:312)))
pbaspect([1 1 1])
colorbar
title('MTF Total','FontSize',14);
xlabel('x-values')
ylabel('y-values')
saveas(gcf,'6-4_MTF Total','png');

figure; set(gca,'FontSize',14);
plot(u(200:312),abs(MTFTotal(257,200:312)))
pbaspect([1 1 1])
title('MTF Total','FontSize',14);
xlabel('u-values')
ylabel('MTF Total')
saveas(gcf,'6-5_MTF Total X','png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% ARTEFACTOS DE IMAGEN %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Planet = imread('Planet.png');
Planet = im2double(Planet);
Planet = rgb2gray(Planet);
NeoPlanet = abs(conv2(PSFTotal,Planet));

Planet = Planet/max(max(Planet));
NeoPlanet = NeoPlanet/max(max(NeoPlanet));

figure; set(gca,'FontSize',14);
subplot(1,2,1);
imagesc(Planet)
pbaspect([1 1 1])
title('Imagen Normal','FontSize',14);
colormap(gray)
subplot(1,2,2); 
imagesc(NeoPlanet(257:770,257:770))
pbaspect([1 1 1])
title('Efecto PSF Total','FontSize',14);
saveas(gcf,'6-6_Efecto en Imagen X','png');

close all