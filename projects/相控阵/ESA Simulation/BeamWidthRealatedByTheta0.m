%% plot beamwidth related by frequency and incident angle
%% input arguements
freq = [1 5 10 15];
lambda = 0.3./freq; 
theta0 = linspace(0, 60, 100); % incident angle
D = 1; % aperture size
K = 0.886; % beamwidth factor

%% calculate beamwidth 
[~, theta0_mat] = meshgrid(lambda, theta0);
bw = lambda*K./(cosd(theta0_mat)*D);
bw_deg = bw*180/pi; % beamwidth in degree;

%% plot
clf
figure(1)
set(gca, 'fontsize', 16);
plot(theta0_mat, bw, 'linewidth', 2);
xlabel('\theta_{0} (degrees)');
ylabel('BW (degress)');
grid off
legend('1GHz', '5GHz', '10GHz', '15GHz');