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
                     
                     
%% RD simulation

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

%% Range Compression
f_rg = linspace(0, 1, SAR.Nrg_cells)*SAR.Fr - SAR.Fr/2;
Ksrc = 2*SAR.Vr^2*SAR.f0^3*(sqrt(1 - SAR.lambda^2*SAR.fdc^2/4/SAR.Vr^2)^3)...
       /(SAR.C*SAR.R0.*SAR.fdc.^2);
Km = SAR.Kr./(1 - SAR.Kr./Ksrc);
Hr_win = createKaiserWin(SAR.Hr_size, SAR.winBeta, SAR.Nrg_cells, round(SAR.Nrg_cells/2));
Hr_match = exp(-1j*pi.*f_rg.^2/Km);
Hr_match = Hr_match.*Hr_win';
Hr_match = ones(SAR.Nrg_lines, 1)*Hr_match;
S_RC = ifft(ifftshift(fftshift(fft(data, [], 2), 2).*Hr_match, 2), [], 2);

%% RCMC
S_RC_RD  = fft(S_RC);
f_az = linspace((SAR.Mamb - 0.5)*SAR.PRF,(SAR.Mamb + 0.5)*SAR.PRF, SAR.Nrg_lines)';
D = (1 -  SAR.lambda^2*f_az.^2/4/SAR.Vr^2).^(1/2);
D = D*ones(1, SAR.Nrg_cells);

t = ones(SAR.Nrg_lines, 1)*(linspace(0, SAR.Nrg_cells - 1, SAR.Nrg_cells))/SAR.Fr;
cell_R0 = (t + SAR.t_start)*SAR.C/2;
deltaR = ((1 - D)./D).*cell_R0;
deltaTau = deltaR*2/SAR.C;

% build interpolation kernals table
interK.oversampleRate = 8;
interK.spacing = 1/interK.oversampleRate;
interK.n = 8;
interK.beta = 2.5;
interK.factor = zeros(interK.oversampleRate, interK.n);
timeSeq = [-4:1:3];


for kk = 1 : interK.oversampleRate
    timeSeq = timeSeq - interK.spacing;
    interK.factor(kk, :) = sinc(timeSeq).*kaiser(interK.n, interK.beta)';
    Norm = sqrt(sum(interK.factor(kk, :).^2));
    interK.factor(kk, :) = interK.factor(kk, :)/Norm;
end

% sinc interpolation
S_RCMC = zeros(SAR.Nrg_lines, SAR.Nrg_cells);
for kaz = 1 : SAR.Nrg_lines
    for krg = 1 : SAR.Nrg_cells
         % quantify deltaR
         tau = t(kaz, krg) + deltaTau(kaz, krg);
         % choose interpolation factors index
         quantK = round(mod(tau, 1/SAR.Fr)/interK.spacing*SAR.Fr);
   
         if quantK == 0
             quantK = interK.oversampleRate;
         end
         
         val = zeros(interK.n, 1);
         ktauL = ceil(tau*SAR.Fr);
         ktauR = ktauL + 1;
         % decrease one to use mod operation and 
         % get cyclical interK.n values of S_RC_RC
         ktauL = ktauL - 1;
         ktauR = ktauR - 1;
         kL = floor(interK.n/2); kR = kL + 1;
         while kL > 0
             ktauL = mod(ktauL + SAR.Nrg_cells, SAR.Nrg_cells);
             val(kL) = S_RC_RD(kaz, ktauL + 1);
             kL = kL - 1;
             ktauL = ktauL - 1;
         end
         while kR <= interK.n
             ktauR = mod(ktauR, SAR.Nrg_cells);
             val(kR) = S_RC_RD(kaz, ktauR + 1);
             kR = kR + 1;
             ktauR = ktauR + 1;
         end
         S_RCMC(kaz, krg) = interK.factor(quantK, :)*val;
    end
end


%% azimuth direction compression
center = find(f_az >= SAR.fdc, 1);
Ha_win = createKaiserWin(SAR.Ha_size, SAR.winBeta, SAR.Nrg_lines,center);
Ha_win = Ha_win*ones(1, SAR.Nrg_cells);
Ha_match = exp(1j*4*pi*cell_R0.*D*SAR.f0/SAR.C); 
% Ha_match = exp(-1j*pi*f_az.^2/SAR.Ka)*ones(1, SAR.Nrg_cells);
Ha_match = Ha_match.*Ha_win;
Sout = ifft(ifftshift(Ha_match.*fftshift(S_RCMC)));
imagesc(abs(Sout));
