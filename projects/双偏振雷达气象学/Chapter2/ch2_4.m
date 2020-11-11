clc;

load('D:/Assignment/Workspace/Radar Assignments/Materials/Chapter2Material/dsddata_20050513.mat');
D = squeeze(dsd_data(:,3,:));    %% Diameter: bin spaced by D 
ND = squeeze(dsd_data(:,6,:));   %% DSD
V = squeeze(dsd_data(:,8,:));     %% velocity
TIME_LEN = size(ND, 2);
D_LEN = size(ND, 1);
times = files(:,8:11);
timeIndex = str2num(times(:,1)).*10+str2num(times(:,2))+(str2num(times(:,3)).*10+str2num(times(:,4)))./60;

% 计算总表面分布、 质量分布、反射率分布
TotalSurfaceD = pi.*D.^2.*ND;
ContentD = (pi/6*10^-6).*D.^3.*ND;
ZD = D.^6.*ND;

% 绘图
figure(1)
subplot(4,1,1)
pcolor(timeIndex(:,1), D(:,1), log10(ND));
colormap(jet)
shading flat
title(strcat('log_{10}(N(D))  /lg(m^{-3}\cdotmm^{-1})'));
ylabel('D  /mm');
colorbar
xlim([7 13])
subplot(4,1,2)
pcolor(timeIndex(:,1), D(:,1), log10(TotalSurfaceD));
colormap(jet)
shading flat
title(strcat('log_{10}(A(D))  /lg(mm^{2}\cdotm^{-3}\cdotmm^{-1})'));
ylabel('D  /mm');
colorbar
xlim([7 13])
subplot(4,1,3)
pcolor(timeIndex(:,1), D(:,1), log10(ContentD));
colormap(jet)
shading flat
title(strcat('log_{10}(M(D))  /lg(g\cdotm^{-3})'));
ylabel('D  /mm');
colorbar
xlim([7 13])
subplot(4,1,4)
pcolor(timeIndex(:,1), D(:,1), log10(ZD));
colormap(jet)
shading flat
title(strcat('log_{10}(Z(D))  /lg(mm^{6}\cdotm^{-3}\cdotmm^{-1})'));
ylabel('D  /mm');
colorbar
xlim([7 13])

% 计算数密度、雨水含量、降雨率、雷达反射率因子
NC = sum(ND, 1);
LWC = sum(ContentD, 1);
R = sum(6*pi*10^-4.*D.^3.*V.*ND, 1);
Z = sum(ZD, 1);   

figure(2)
subplot(4,1,1)
plot(timeIndex(:,1), NC);
title(strcat('Number Concentration  /m^{-3}'));
xlim([7 13])
subplot(4,1,2)
plot(timeIndex(:,1), LWC);
title(strcat('Liquid Water Content  /g\cdotm^{-3}'));
xlim([7 13])
subplot(4,1,3)
plot(timeIndex(:,1), R);
title(strcat('Rain Rate  /mm\cdoth^{-1}'));
xlim([7 13])
subplot(4,1,4)
plot(timeIndex(:,1), Z);
title(strcat('Reflectivity  /mm^{6}\cdotm^{-3}\cdotmm^{-1}'));
xlim([7 13])

M2 = sum(D.^2.*ND, 1);
M3 = sum(D.^3.*ND, 1);
M4 = sum(D.^4.*ND, 1);
M6 = sum(D.^6.*ND, 1);
%M-P Model
MP_N0 = 8*10^6.*ones(size(M2));
MP_lambda = (gamma(4)*MP_N0./M3).^(1/4);
MP_Model = MP_N0.*exp(-MP_lambda.*D);
%EXP Model
EXP_lambda = ((M2.*gamma(5))./(M4.*gamma(3))).^(1/2);
EXP_N0 = M2.*(EXP_lambda.^3)/gamma(3);
EXP_Model = EXP_N0.*exp(-EXP_lambda.*D);
%GAMMA Model
GAMMA_eta = M4.^2./M2./M6;
GAMMA_miu = ((7 - 11*GAMMA_eta) - (GAMMA_eta.^2 + 14*GAMMA_eta + 1).^(1/2))./(2.*(GAMMA_eta - 1));
GAMMA_lambda = (M2./M4.*(GAMMA_eta + 3).*(GAMMA_eta + 4)).^(1/2);
GAMMA_N0 = M2.*(GAMMA_lambda.^(GAMMA_miu + 3))./gamma(GAMMA_miu + 3);
GAMMA_Model = GAMMA_N0.*(D.^GAMMA_miu).*exp(-GAMMA_lambda.*D);

% 绘图
figure(3)
subplot(3,1,1)
pcolor(timeIndex(:,1), D(:,1), log10(EXP_Model));
colormap(jet)
shading flat
title(strcat('EXP log_{10}(N(D))  /lg(m^{-3}\cdotmm^{-1})'));
ylabel('D  /mm');
colorbar
xlim([7 13])
subplot(3,1,2)
pcolor(timeIndex(:,1), D(:,1), log10(MP_Model));
colormap(jet)
shading flat
title(strcat('M-P log_{10}(N(D))  /lg(m^{-3}\cdotmm^{-1})'));
ylabel('D  /mm');
colorbar
xlim([7 13])
subplot(3,1,3)
pcolor(timeIndex(:,1), D(:,1), log10(GAMMA_Model));
colormap(jet)
shading flat
title(strcat('GAMMA log_{10}(N(D))  /lg(m^{-3}\cdotmm^{-1})'));
ylabel('D  /mm');
colorbar
xlim([7 13])

% % 第35个滴谱拟合模型曲线
figure(4)
subplot(3,1,1)
plot(D(:,35), MP_Model(:, 35));
title(strcat('M-P log_{10}(N(D))  /m^{-3}\cdotmm^{-1}'));
subplot(3,1,2)
plot(D(:,35), EXP_Model(:, 35));
title(strcat('EXP log_{10}(N(D))  /m^{-3}\cdotmm^{-1}'));
subplot(3,1,3)
plot(D(:,35), GAMMA_Model(:, 35));
title(strcat('GAMMA log_{10}(N(D))  /m^{-3}\cdotmm^{-1}'));

% 利用计算出来的模型来估计数密度、含水量、降雨率、反射率因子
calculateNC_W_R_Z(timeIndex, D, MP_Model, V, 'M-P');
calculateNC_W_R_Z(timeIndex, D, EXP_Model, V, 'EXP');
calculateNC_W_R_Z(timeIndex, D, GAMMA_Model, V, 'GAMMA');