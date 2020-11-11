function computeAzSpectra(data, SAR)
%% compute azimuth spectra
set( 0, 'DefaultTextFontSize',   12 )  % Plotting defaults
set( 0, 'DefaultLineLineWidth', 1.5 )
set( 0, 'DefaultAxesFontSize',    8 )
lfs = 11;

disp ' '
disp '---------------------------------------------------------'
fprintf(' UBC RRSG - Plot the azimuth spectrum of each data block')
disp ' '
disp '---------------------------------------------------------'

Nrowsg = 3;     % Number of subplots in row direction of the figure 
Ncolsg = 3;     % Number of subplots in column direction of the figure 
Nspect = Nrowsg*Ncolsg;         % Total number of spectra to calculate
Nrglpb = floor( SAR.Nrg_cells/Nspect );  %  No. of range cells per spectra
wd = 0.81/Ncolsg;   dx = wd + 0.045;   x0 = 0.07;    
ht = 0.39/Nrowsg;   dy = 0.28;         y0 = 1-dy;  

%% test azimuth spectra
disp 'Compute and plot azimuth power spectra'
tic
figure(202), clf
freq = [0:SAR.Nrg_lines-1]*SAR.PRF/SAR.Nrg_lines-1;

%  Find azimuth power spectra and average
ysc = zeros(Nspect, 1);
DATA_aver = zeros(Nspect, SAR.Nrg_lines);
for krg = 1 : Nspect
    r1 = 1 + (krg-1)*Nrglpb;   r2 = r1 + Nrglpb - 1;
    DATA = fft( data(:,r1:r2) ); % azimuth direction fft
    DATA_aver(krg,:) = mean( abs(DATA.').^2 )/1000000;
    ysc(krg) = 1.05*max(DATA_aver(krg,:));
end  % of for krg = 1 : Nspect
ysc0 = max( ysc );      %  Common vertical scaling for all the spectra

Ffrac = zeros(Nspect, 1);
for krg = 1 : Nspect
    subplot(Nrowsg, Ncolsg, krg)
    plot( freq, DATA_aver(krg,:) ),   grid,   hold on
    set( gca, 'Pos',...
       [x0+dx*mod((krg-1),Ncolsg)  y0-dy*floor((krg-1)/Ncolsg)  wd ht])
    axis([0 SAR.PRF  0 ysc0]);

    azim_spec = fft( DATA_aver(krg,:) )/ SAR.Nrg_lines;% azimuth 
    angle_first_harmonic = -angle( azim_spec(2) );
    Ffrac(krg) = SAR.PRF * angle_first_harmonic / (2*pi);
    if Ffrac(krg) < 0,   Ffrac(krg) = Ffrac(krg) + SAR.PRF;   end
    sine_fit = real(azim_spec(2)) * cos(2*pi*freq/SAR.PRF) - ...
            imag(azim_spec(2)) * sin(2*pi*freq/SAR.PRF) + 0.5*azim_spec(1);
    plot( freq, 2*sine_fit, 'r--' )
    
    if krg > Nspect - Ncolsg
        xlabel('Azimuth frequency  (Hz)  \rightarrow', 'FontSize', lfs )
    end
    if mod(krg,Ncolsg) == 1
        ylabel('Power  \rightarrow', 'FontSize', lfs )
    end
end
end
