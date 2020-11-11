%% RADARSAT_RDA_Imaging.m
% ��RD�㷨��RADARSAT-1�Ļز����ݳ���
% author: ��Ժ��� ����(07042185)
% Ϊ�˽�ʡ������ڴ濪���������а����˽϶��clear���
% Ϊ�˼��ٳ�������ʱ���ڴ濪����ֵ�������ʽ�ֶಽʵ��

%% ��ȡ����
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
                     
                     
% ���ļ�specify_parameters.m��ָ����ز���
specify_parameters;
% ���ļ��г�ȡ����
extract_data;
% �����ݿ�1���г���
AzimuthBlock=1;
% ��������
load(strcat(output_path,'CDdata',num2str(AzimuthBlock)));
file_pre = strcat( output_path, output_prefix, '_', num2str(AzimuthBlock) );
% �õ��Զ��������(AGC)ֵ
AGC_values = load_AGC_block( file_pre, first_rg_line, ...
    Nrg_lines_blk, AzimuthBlock , UseMATfiles );
% AGC��ԭ���õ���ʵ�Ļز�����
org_data = load_DATA_block( file_pre, output_path, Nrg_lines_blk, ...
    Nrg_cells, AGC_values, AzimuthBlock, UseMATfiles );
% �������е�clear��䶼��Ϊ�˽�ʡ�ڴ濪��
clear AGC_values

%% ��һ����ȷRADARSAT-1����ز���

length_replica  =  2880;         % Total length (I&Q) of replica record
tot_Nrg_cells   =  9288;         % Total number of range cells per line
tot_Nrg_lines   = 19432;         % Total number of range lines (records)
first_replica   =     7;         % First record that contains the replica
PRF             = 1256.98;       % Pulse Reputation Frequency (Hz)
Fa              = PRF;           % ��λ�����Ƶ��
Fr              = 32.317e+6;     % Radar sampling rate (Hz)
f0              = 5.300e+9;      % Radar center frequency (Hz)
c               = 2.9979e+8;     % Speed of light (m/s)
t0              = 0.0065956;     % data window start time
R0              = 0.0065956*c/2; % Slant range of first radar sample (m)
Nrepl           = 1349;          % No. of valid samples in the replica
Kr              = 0.72135e+12;   % FM rate of radar pulse (Hz/s)
Tr              = 41.75e-6;      % Chirp duration (s)
Vr              = 7062;          % Effective radar velocity (m/s)
WaveLength      = c/f0;          % Radar wavelength (m)
Ka              = 1733;          % Azimuth FM rate (Hz/s)
f_dc            = -6900;         % Doppler centroid (Hz)
Lr              = 15;            % Range Aperture
La              = 1.5;           % Azimuth aperture
Na_mf           = 705;           % ��λ����ƥ���˲�������
Ta              = Na_mf/Fa;      % the duration of the azimuth match filter

%% ��RADARSAT-1�Ĳ�����һ���õ��������

% �õ�������Ļز�����
data=(org_data); %(1:2048,:);)
clear org_data;
% NaΪrange line����NrΪrange��������
[Na, Nr]=size(data); 

% �����������Ƶ�ʷֱ���
d_fr=Fr/Nr;
% ���㷽λ����Ƶ�ʷֱ���
d_fa=Fa/Na;

% �����������ʱ�����
RangeTimeLimit=Nr/Fr;
RangeTime=(0:Nr-1)/Fr;

% ���㷽λ����ʱ�����
AzimuthTimeLimit=Na/Fa;
AzimuthTime=(0:Na-1)/Fa;

%% (1)RC��RCMC
data=fft(data,[],2); % �������ϵ�FFT

% �����������ϵ�����ź�
% �������ϵ�LFM�ź�
s_r=exp(-i*pi*Kr*(RangeTime-Tr/2).^2).*(abs(RangeTime-Tr/2)<=Tr/2);
% ��������LFM�źŵ�Ƶ��
s_r=fft(s_r);
% ��������ƥ���˲���Ƶ��
s_r=conj(s_r);
s_r=ones(Na,1)*s_r; % ����

data=data.*s_r; % RC
clear s_r
% ����RCMC�ź�
fr=d_fr*(0:Nr-1); % �������Ƶ�ʲ����� (Hz)
s_rcmc=exp(-i*4*pi*f_dc*WaveLength/2*( (AzimuthTime)'*ones(1,Nr) )/c.*( ones(Na,1)*fr ));
clear fr
data=data.*s_rcmc; % RCMC
clear s_rcmc

% IFFTתΪ�������ϵ�ʱ���ź�
data=ifft(data,[],2);
%% (2) AC
data=fft(data,[],1); % ��λ���ϵ�FFT

% ������λ���ϵ�����ź�
% ��λ��Ƶб�����������ԱȽ�ǿʱ��uncomment��������
R=c*RangeTime/2+c*first_rg_cell/Fr/2+R0; % ��ʱ���Ӧ��б��
K=-2*Vr^2/WaveLength./R; % б���Ӧ�ĵ�Ƶб��
clear R
clear RangeTime AzimuthTime

% ��λ��Ƶб�����������ԱȽ���ʱ���ù̶���Ka����AC
K=-Ka*ones(1,Nr);

% ���㷽λ����ƥ���˲�����Ӧ��ʱ��
AzimuthMatchFilterTime=(0:Na_mf-1)/Fa; 
% ��λ���ϵ�LMF�ź�
s_a=exp(-i*2*pi*f_dc*((AzimuthMatchFilterTime-Ta/2)'*ones(1,Nr))+...
    i*pi*((AzimuthMatchFilterTime-Ta/2)'.^2)*K).*...
    (abs((AzimuthMatchFilterTime-Ta/2)'*ones(1,Nr))<=Ta/2); 
clear K AzimuthMatchFilterTime
% ������Na�㳤
s_a=[s_a;zeros(Na-Na_mf,Nr)]; 
% ��λ����LFM�źŵ�Ƶ��
s_a=fft(s_a,[],1); 
% ��λ����ƥ���˲�����Ƶ��
s_a=conj(s_a); 

data=data.*s_a; % AC
clear s_a;

% IFFTתΪ��λ���ϵ�ʱ���ź�
data=ifft(data,[],1);


%%
% ======================== The End ========================