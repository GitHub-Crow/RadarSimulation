%% compute subarray architecture pattern

%% init paramenters
IBWratio = 1.1;
f0 = 4;
f = IBWratio*f0;
lambda = 0.3/f;
lambda0 = 0.3/f0;
d = lambda/2;

theta = linspace(-10, 50, 721); theta0 = 30;
u = sind(theta); u0 = sind(theta0);

OSA.ratio = 2;
EF = 1.5;

SA.nele = 8; SA.nsas = 20;
OSA.nele = OSA.ratio*SA.nele; % number of elements in OSA
OSA.nsas = SA.nsas - (OSA.ratio - 1); % number of OSAs in backend OSA AF
%% init weights
Taylor.nbar = 6;
Taylor.sll = 35;
Taylor.set = true;
SA.wgts.ele = ones(1, SA.nele);
SA.wgts.sas = ones(1, SA.nsas);
OSA.wgts.ele = ones(1, OSA.nele);
OSA.wgts.sas = ones(1, OSA.nsas);

if Taylor.set
    OSA.wgts.ele = getTaylorWeights(OSA.nele, Taylor.sll, Taylor.nbar);
    SA.wgts.ele = getTaylorWeights(SA.nele, Taylor.sll, Taylor.nbar);
end
%% compute pattern
EP = compute1D_EP(theta, EF);
SA.AFe = compute1D_AF(SA.wgts.ele, SA.nele, d, f0, f, u0, u);
OSA.AFe = compute1D_AF(OSA.wgts.ele, OSA.nele, d, f0, f, u0, u);
% backend beamforming
SA.AF = compute1D_AF(SA.wgts.sas, SA.nsas, d, f, f, u0, u);
OSA.AF = compute1D_AF(OSA.wgts.sas, OSA.nsas, d, f, f, u0, u);

SA.PAT = EP.*SA.AFe.*SA.AF;
OSA.PAT = EP.*OSA.AFe.*OSA.AF;
[~, ~, EP_dBNorm] = processVector(EP);
[~, ~, SA.AF_dBNorm] = processVector(SA.AF);
[~, ~, SA.AFe_dBNorm] = processVector(SA.AFe);
[~, ~, OSA.AF_dBNorm] = processVector(OSA.AF);
[~, ~, OSA.AFe_dBNorm] = processVector(OSA.AFe);
[~, ~, SA.PAT_dBNorm] = processVector(SA.PAT);
[~, ~, OSA.PAT_dBNorm] = processVector(OSA.PAT);
% for normalized back end pattern rather than output pattern.
% compensate Array.PAT value loss caused by log operator 
OSA.compensation = EP_dBNorm(u == u0) + OSA.PAT_dBNorm(u == u0);
SA.compensation = EP_dBNorm(u == u0) + SA.PAT_dBNorm(u == u0);
%% plot
figure(1); clf
set(gcf,'DefaultLineLineWidth',2)
plot(theta, OSA.AFe_dBNorm, 'Color', 'b');
hold on
plot(theta, OSA.AF_dBNorm, 'Color', 'k');
hold on
plot(theta, OSA.PAT_dBNorm + OSA.compensation, 'Color', 'g');
hold on
set(gca, 'fontsize', 14, 'fontweight', 'bold');
title(['\it Composite Array Pattern with ',num2str(OSA.ratio),...
': 1 Overlap and f = ',num2str(IBWratio),'*f_o'])
xlabel('\theta (degree)'); ylabel('dB');
legend('OSA Subarray Pattern', 'OSA Backend AF', 'OSA pattern');
axis ([theta(1), theta(end), -80, 0]);

figure(2); clf
set(gcf,'DefaultLineLineWidth',2)
plot(theta, SA.AFe_dBNorm, 'Color', 'b');
hold on
plot(theta, SA.AF_dBNorm, 'Color', 'k');
hold on
plot(theta, SA.PAT_dBNorm + SA.compensation, 'Color', 'g');
hold on
set(gca, 'fontsize', 14, 'fontweight', 'bold');
title(['\it Composite Array Pattern with No Overlap and f = ',num2str(IBWratio),'*f_o'])
xlabel('\theta (degree)'); ylabel('dB');
legend('SA Subarray Pattern', 'SA Backend AF', 'SA pattern');
axis ([theta(1), theta(end), -80, 0]);
