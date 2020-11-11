%% 1D pattern
close all; clear all;
%% input parameters
Array.f = 10; % operating frequency in GHz
Array.f0 = 10; % tune frequency in GHz
Array.nele = 30; % number of array element
Array.d = 0.5*(0.3/Array.f0); % array element distance
Array.EF = 1.35;
Array.weightFlag = 1; % 0 = uniform, 1 = taylor weighting
% Taylor Weighting paramenters
Array.Taylor.Nbar = 5;
Array.Taylor.SLL = 30;

% theta angle paramenters
theta.scan = 0;
theta.minAngle = -90;
theta.maxAngle = 90;
theta.nAngles = 721;

%% set pattern of amplitude weight
if Array.weightFlag == 0
    Array.amWeight = ones(1, Array.nele);
else
    Array.amWeight = getTaylorWeights(Array.nele, Array.Taylor.SLL, Array.Taylor.Nbar);
end

theta.vec = linspace(theta.minAngle, theta.maxAngle, theta.nAngles);
theta.u = sind(theta.vec);
theta.u0 = sind(theta.scan);

% init element pattern, array factor, pattern
Array.size = size(theta.vec);
Array.EP = zeros(Array.size);
Array.AF = zeros(Array.size);
Array.PAT = zeros(Array.size);

% compute element pattern, array factor, pattern
Array.EP = compute1D_EP(theta.vec, Array.EF);
Array.AF = compute1D_AF(Array.amWeight, Array.nele, Array.d, ...
           Array.f0, Array.f, theta.u0, theta.u);
Array.PAT = compute1D_PAT(Array.EP, Array.AF);

%% plot
[EP.mag, EP.mag_dB, EP.mag_dB_Norm] = processVector(Array.EP);
[AF.mag, AF.mag_dB, AF.mag_dB_Norm] = processVector(Array.AF);
[PAT.mag, PAT.mag_dB, PAT.mag_dB_Norm] = processVector(Array.PAT);

clf
figure(1)
set(gca, 'fontsize', 14, 'fontweight', 'bold');

plot(theta.vec, EP.mag_dB_Norm, 'Color', 'r', 'linestyle', ':', 'linewidth', 2);
hold on
plot(theta.vec, AF.mag_dB_Norm, 'Color', 'b', 'linestyle', ':',  'linewidth', 2);
hold on
plot(theta.vec, PAT.mag_dB_Norm, 'Color', 'k', 'linestyle', '-',  'linewidth', 3);
hold on
xlabel('\theta (degree)'); ylabel('dB');
legend('EP', 'AF', 'PAT');
axis ([theta.minAngle, theta.maxAngle, -120 20]);