clc; close all;
fontsize = 16;

%% Main : clear ground phase
phase = imread('pha.tif');
TD = exp(1j.*phase);
[Na, Nr] = size(TD);

%% clear azimuth stripe frequency
FDA = fft(TD, [], 2);
[~, index] = max(sum(FDA, 1));                                % find stripe frequency
FDA_new = horzcat(FDA(:, index : Nr), FDA(:, 1 : index - 1)); % move to zero frequency
phase_azimuthGC = angle(ifft(FDA_new, [], 2));                % get phase image

%% rebuild time sequence
TD__azimuthGC = exp(1j*phase_azimuthGC);

%% clear range stripe frequency
FDR = fft(TD__azimuthGC);  % range fft
[~, index] = max(sum(FDR, 2));                                % find stripe frequency
FDR_new = vertcat(FDR(index : Na, :), FDR(1 : index - 1, :)); % move to zero frequency
phase_GC = angle(ifft(FDR_new));                              % get phase image

%% plot original phase image and processed phase image
figure ;
subplot(121)
imshow(phase)
title('original phase image', 'fontsize', fontsize);
subplot(122)
imshow(phase_GC);
title('pahse image cleared ground phase', 'fontsize', fontsize);
imwrite(phase_GC, 'phase_GC.png');

%% Main : location

%% init paramenters
R0_main = 1443640.2;                            % main image 
PS_main = [1466309.6  -4601880.8  4914932.0];   % main image location
VS_main = [-61.062392  -5620.2862  -5236.1316];

R0_assist = 1443566.7;
PS_assist = [1466545.1  -4601144.2  4915020.0];

lambda = 0.031250000;
fdc = -295.19346;
PT0 = [884704.82  -4306714.2  4605148.9];

%% solve R-D equations
syms PT [1 3]
eqn1 = R0_main/2 == sqrt(sum((PT - PS_main).^2));
eqn2 = R0_assist - R0_main/2 == sqrt(sum((PT - PS_assist).^2));
eqn3 = fdc == -2*(sum(VS_main.*(PS_main - PT)))/lambda/R0_main;

eqns = [eqn1 eqn2 eqn3];
ANS = solve(eqns, PT);
ANS = struct2cell(ANS);

%% choose the closest answer from PT0 
PT = zeros(1, 3);
curPT = zeros(1, 3);
ClosestDistFromPT0 = intmax('int64');

N = length(ANS{1});
for kk = 1 : N
    for ii = 1 : 3
        curPT(ii) = eval(ANS{ii}(kk));
    end
    dist = sum((curPT - PT0).^2);
    if dist < ClosestDistFromPT0
        PT = curPT;
        ClosestDistFromPT0 = dist;
    end
end
log = sprintf('PT located at [%.2f  %.2f  %.2f].', PT(1), PT(2), PT(3));

disp(log);