%% compute subarray architecture pattern

%% init parameters
IBWratio = 1.1;
f0 = 4;
f = IBWratio*f0;
lambda = 0.3/f;
lambda0 = 0.3/f0;
d = lambda/2;

theta = linspace(-90, 90, 721); theta0 = 30;
u = sind(theta); u0 = sind(theta0);
Subarray.nele = 4;
Array.nele = 4;

EF = 1.5;
%% init weights
Subarray.wgts = ones(1, Subarray.nele);
Array.wgts = ones(1, Array.nele);

%% compute pattern
Subarray.AF = compute1D_AF(Subarray.wgts, Subarray.nele, d, f0, f, u0, u);
% backend array
Array.EP = compute1D_EP(theta, EF);
Array.AF = compute1D_AF(Array.wgts, Array.nele, d, f, f, u0, u);
Array.PAT = Array.EP.*Subarray.AF.*Array.AF;

[~, ~, Subarray.AF_dBNorm] = processVector(Subarray.AF);
[~, ~, Array.AF_dBNorm] = processVector(Array.AF);
[~, ~, Array.EP_dBNorm] = processVector(Array.EP);
[~, ~, Array.PAT_dBNorm] = processVector(Array.PAT);
% for normalized back end pattern rather than output pattern.
% compensate Array.PAT value loss caused by log operator 
compensation = Array.EP_dBNorm(u == u0) + Subarray.AF_dBNorm(u == u0);
%% plot
clf
figure(1)
set(gcf,'DefaultLineLineWidth',2)
plot(theta, Array.EP_dBNorm, 'Color', 'r');
hold on
plot(theta, Subarray.AF_dBNorm, 'Color', 'b');
hold on
plot(theta, Array.AF_dBNorm, 'Color', 'k');
hold on
plot(theta, Array.PAT_dBNorm + compensation, 'Color', 'g');
hold on

set(gca, 'fontsize', 14, 'fontweight', 'bold');
xlabel('\theta (degree)'); ylabel('dB');
legend('EP', 'Subarray.AF', 'AF', 'PAT');
axis ([theta(1), theta(end), -80, 0]);
