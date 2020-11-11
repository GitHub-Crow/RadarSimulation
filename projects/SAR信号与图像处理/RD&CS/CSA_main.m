%% load SAR data from CD
clear,    format compact
set( 0, 'DefaultTextFontSize',   12 )  % Plotting defaults
set( 0, 'DefaultLineLineWidth', 1.5 )
set( 0, 'DefaultAxesFontSize',    8 )
load CD_run_params
block = 1;  
file_pre = strcat( output_path, output_prefix, '_', num2str(block) );

disp ' '
disp (['Load or Extract AGC setting and Data for block ' num2str(block) ])
%  Load a block of 'AGC_values'
AGC_values = load_AGC_block( file_pre, first_rg_line, ...
                                      Nrg_lines_blk, block , UseMATfiles );

%  Load a block of raw SAR data
data = load_DATA_block( file_pre, output_path, Nrg_lines_blk, ...
                         Nrg_cells, AGC_values, block, UseMATfiles );
                     
                     
%% CS simulation

%% init SAR paramenters
SAR.f0 = f0;
SAR.lambda = c/SAR.f0;
SAR.Kr = Kr;
SAR.R0 = R0;
SAR.Nrg_cells = Nrg_cells;
SAR.Nrg_lines = Nrg_lines_blk;
SAR.Vr = 7062;
SAR.PRF = PRF;
SAR.Fr = Fr;
SAR.Tr = Tr;
SAR.C = 299790000;
SAR.t_start = 0.0065956;
SAR.fdc = -6900;
SAR.Ka = 1733;
SAR.Mamb = round(SAR.fdc/SAR.PRF);
SAR.winBeta = 2.5;
SAR.Ha_size = 705;
SAR.Hr_size = round(SAR.Tr*SAR.Fr);

%% chrip scaling
t = ones(SAR.Nrg_lines, 1)*(linspace(0, SAR.Nrg_cells - 1, SAR.Nrg_cells))/SAR.Fr;
cell_R0 = (t + SAR.t_start)*SAR.C/2;

f_az = linspace((SAR.Mamb - 0.5)*SAR.PRF,(SAR.Mamb + 0.5)*SAR.PRF, SAR.Nrg_lines)';
f_az = f_az*ones(1, SAR.Nrg_cells);
D = (1 -  SAR.C^2*f_az.^2/(4*SAR.Vr^2*SAR.f0^2)).^(1/2);
Km = SAR.Kr./(1 - SAR.Kr*SAR.C*cell_R0.*f_az.^2./(2*SAR.Vr^2*SAR.f0^3*D.^3));
f_ref = SAR.fdc;
Dref = (1 -  SAR.C^2*f_ref^2/(4*SAR.Vr^2*SAR.f0^2))^(1/2);
Rref = (SAR.t_start + SAR.Nrg_cells/2/SAR.Fr)*SAR.C/2;
t_re = 2/SAR.C*(cell_R0./Dref - Rref/Dref);
S_sc = exp(1j*pi*Km.*(Dref./D - 1).*t_re.^2);
S = ifftshift(fftshift(fft(data)).*S_sc);

%% Range Compression
S_2DF = fftshift(fft(S, [], 2), 2);
f_rg = linspace(0, 1, SAR.Nrg_cells)*SAR.Fr - SAR.Fr/2;
f_rg = ones(SAR.Nrg_lines, 1)*f_rg;

Hr_match = exp(-1j*pi*D./Km./Dref.*f_rg.^2).*exp(1j*4*pi/SAR.C*(1./D - 1/Dref).*Rref.*f_rg);
S_2DF = ifftshift(S_2DF.*Hr_match, 2);

%% Azimuth Conpression
S_RD = ifftshift(ifft(S_2DF, [], 2));
Ha_match = exp(1j*4*pi*cell_R0.*D*SAR.f0/SAR.C).*...
           exp(1j*4*pi*Km/SAR.C^2.*(1 - D/Dref).*(cell_R0./D - Rref./D).^2); 
center = find(f_az >= SAR.fdc, 1);
Ha_win = createKaiserWin(SAR.Ha_size, SAR.winBeta, SAR.Nrg_lines,center);
Ha_win = Ha_win*ones(1, SAR.Nrg_cells);
Ha_match = Ha_match.*Ha_win;
Sout = ifft(ifftshift(Ha_match.*S_RD));
colormap jet
imagesc(abs(Sout));
