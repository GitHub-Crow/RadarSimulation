%% Code to compute subarrayed architecture pattern
% Arik D. Brown
%% Input Parameters
IBWratio=1.1;%IBWratio -> f=fo
fo=4;%GHz Tune Frequency
f=IBWratio*fo;%GHz Operating Frequency
lambda=0.3/f;%inches
lambdao=0.3/fo;%inches
d=lambdao/2;%inches
theta=linspace(-90,90,721);%deg
thetao=30;%deg98 Electronically Scanned Arrays: MATLAB? Modeling and Simulation
u=sind(theta);
uo=sind(thetao);
SA.nelems=4;%Number of elements in Subarray
AF.nelems=4;%Number of elements in Backend AF
EF=1.5;
SA.wgts=ones(1,SA.nelems);
AF.wgts=ones(1,AF.nelems);
plotfigs.flags.EP=1;
plotfigs.flags.SA=1;
plotfigs.flags.AF=1;
plotfigs.flags.PAT=1;
plotfigs.flags.ALL=1;
plotfigs.axis.xlims=[-90 90];
plotfigs.axis.ylims=[-40 0];
%% Compute Pattern
% Element Pattern
EP = compute1D_EP(theta,EF);
[EP_mag, EP_dB, EP_dBnorm] = processVector(EP);
% Subarray AF
%!!! To simulate an array of elements without phase shifters

%the last input to Compute_1D_AF
SA.AF = compute1D_AF(SA.wgts,SA.nelems,d,fo,f,uo,u);
[SA.AF_mag, SA.AF_dB, SA.AF_dBnorm] = processVector(SA.AF);
%Backend AF
AF.AF = compute1D_AF(AF.wgts,AF.nelems,d,fo,f,uo,u);
[AF.AF_mag, AF.AF_dB, AF.AF_dBnorm] = processVector(AF.AF);
%Pattern = Element Pattern x Subarray AF Pattern x AF

PAT=EP.*SA.AF.*AF.AF;
[PAT_mag,PAT_dB,PAT_dBnorm] = processVector(PAT);
SA.scanvalue.afsa=SA.AF_dBnorm(u==uo);
EPscanvalue=EP_dBnorm(u==uo);
patnorm.sa=SA.scanvalue.afsa+EPscanvalue;

%Plot Pattern in dB, Unnormalized
figure,clf
set(gcf,'DefaultLineLineWidth',1.5)
set(gcf,'DefaultTextFontSize',12,'DefaultTextFontWeight','bold')
plot(theta,EP_dBnorm,'color',[0 0 0]),hold
plot(theta,SA.AF_dBnorm,'color',[0 .7 0])
plot(theta,AF.AF_dBnorm,'color',[.7 0 1])
plot(theta,PAT_dBnorm+patnorm.sa,'color',[0 0 1])
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
title(['Composite Array Pattern, f =',num2str(IBWratio),'*f_{o}'],...
'FontSize',14,'FontWeight','bold')
xlabel('\theta (degrees)','FontSize',12,'FontWeight','bold')
ylabel('dB','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
set(gcf,'color','white')
legend('EP','SA.AF','AF','Pattern=EP*SA.AF*AF')