close all;
%% first order phase compensation
Rref = near_range + Nr/2/Fr*C/2; 
alpha = acos(ref_height(1)/Rref);                
deltaR = (height - ref_height)*cos(alpha) + cross*sin(alpha);  
deltaR_Mat = deltaR'*ones(1, Nr);
S = data.*exp(1j*4*pi*deltaR_Mat/lambda);

%% range direction resample
fr = linspace(0, Fr, Nr);
fr_Mat = ones(Na, 1)*fr;
phase_a = exp(1j*4*pi*fr_Mat.*deltaR_Mat/C);
S = ifft(fft(S.*phase_a, [], 2), [], 2);

%% azimuth direction resample
ta = linspace(0, Na - 1, Na)/PRF;
t_moco = time - time(1);
% azimuth direction resample is omited cause
% interval of samples is identical to MOCO

%% chrip scaling
t = ones(Na, 1)*(linspace(0, Nr - 1, Nr))/Fr;
cell_R0 = near_range + t*C/2;
fdc = 2*Vr*sin_theta/lambda;
Mamb = round(fdc/PRF);
fa = linspace((Mamb - 0.5)*PRF,(Mamb + 0.5)*PRF, Na)';
fa_Mat = fa*ones(1, Nr);
D = (1 -  C^2*fa_Mat.^2/(4*Vr^2*f0^2)).^(1/2);
Km = Kr./(1 - Kr*C*cell_R0.*fa_Mat.^2./(2*Vr^2*f0^3*D.^3));
fref = fdc;
Dref = (1 -  C^2*fref^2/(4*Vr^2*f0^2))^(1/2);

t_re = 2/C*(cell_R0./Dref - Rref/Dref);
S_cs = exp(1j*pi*Km.*(Dref./D - 1).*t_re.^2);
S_RD = ifftshift(fftshift(fft(data)).*S_cs);

%% Range Compression
S_2DF = fftshift(fft(S_RD, [], 2), 2);
Hr_match = exp(-1j*pi*D./Km./Dref.*fr_Mat.^2).*exp(1j*4*pi/C*(1./D - 1/Dref).*Rref.*fr_Mat);
S_2DF = ifftshift(S_2DF.*Hr_match, 2);
S_2DT = ifft(ifft(S_2DF, [], 2));

%% second order phase compensation
cross_Mat = cross'*ones(1, Nr);
height_Mat = height'*ones(1, Nr);
cell_R0_real = sqrt(((cell_R0.^2 - ref_height(1)^2).^(1/2) - cross_Mat).^2 + height_Mat.^2);
Rref_real = sqrt(((Rref.^2 - ref_height(1)^2).^(1/2) - cross_Mat).^2 + height_Mat.^2);
dR = (cell_R0_real - cell_R0) - (Rref_real - Rref);
S_2DT = S_2DT.*exp(1j*4*pi*dR/lambda);

%% Azimuth Conpression
S_RD = fftshift(ifft(S_2DT));

Ha_size = Na/2; winBeta = 2.5;
Ha_match = exp(1j*4*pi*cell_R0.*D*f0/C).*...
           exp(1j*4*pi*Km/C^2.*(1 - D/Dref).*(cell_R0./D - Rref./D).^2); 
center = find(fa >= fdc, 1);
Ha_win = createKaiserWin(Ha_size, winBeta, Na, center);
Ha_win = Ha_win*ones(1, Nr);
Ha_match = Ha_match.*Ha_win;
Sout = ifft(ifftshift(Ha_match.*S_RD));
colormap jet
imagesc(abs(Sout));
