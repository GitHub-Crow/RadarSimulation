%% RADARSAT_RDA_Imaging.m
% 用RD算法对RADARSAT-1的回波数据成像
% author: 四院五队 王建(07042185)
% 为了节省程序的内存开销，代码中包含了较多的clear语句
% 为了减少程序运行时的内存开销峰值，长表达式分多步实现

%% 抽取数据
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
                     
                     
% 在文件specify_parameters.m中指定相关参数
specify_parameters;
% 从文件中抽取数据
extract_data;
% 对数据块1进行成像
AzimuthBlock=1;
% 加载数据
load(strcat(output_path,'CDdata',num2str(AzimuthBlock)));
file_pre = strcat( output_path, output_prefix, '_', num2str(AzimuthBlock) );
% 得到自动增益控制(AGC)值
AGC_values = load_AGC_block( file_pre, first_rg_line, ...
    Nrg_lines_blk, AzimuthBlock , UseMATfiles );
% AGC还原，得到真实的回波数据
org_data = load_DATA_block( file_pre, output_path, Nrg_lines_blk, ...
    Nrg_cells, AGC_values, AzimuthBlock, UseMATfiles );
% 以下所有的clear语句都是为了节省内存开销
clear AGC_values

%% 进一步明确RADARSAT-1的相关参数

length_replica  =  2880;         % Total length (I&Q) of replica record
tot_Nrg_cells   =  9288;         % Total number of range cells per line
tot_Nrg_lines   = 19432;         % Total number of range lines (records)
first_replica   =     7;         % First record that contains the replica
PRF             = 1256.98;       % Pulse Reputation Frequency (Hz)
Fa              = PRF;           % 方位向采样频率
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
Na_mf           = 705;           % 方位向上匹配滤波器长度
Ta              = Na_mf/Fa;      % the duration of the azimuth match filter

%% 由RADARSAT-1的参数进一步得到所需参数

% 得到待成像的回波数据
data=(org_data); %(1:2048,:);)
clear org_data;
% Na为range line数，Nr为range采样点数
[Na, Nr]=size(data); 

% 计算距离向上频率分辩率
d_fr=Fr/Nr;
% 计算方位向上频率分辩率
d_fa=Fa/Na;

% 计算距离向上时间变量
RangeTimeLimit=Nr/Fr;
RangeTime=(0:Nr-1)/Fr;

% 计算方位向上时间变量
AzimuthTimeLimit=Na/Fa;
AzimuthTime=(0:Na-1)/Fa;

%% (1)RC和RCMC
data=fft(data,[],2); % 距离向上的FFT

% 产生距离向上的相关信号
% 距离向上的LFM信号
s_r=exp(-i*pi*Kr*(RangeTime-Tr/2).^2).*(abs(RangeTime-Tr/2)<=Tr/2);
% 距离向上LFM信号的频谱
s_r=fft(s_r);
% 距离向上匹配滤波器频谱
s_r=conj(s_r);
s_r=ones(Na,1)*s_r; % 矩阵化

data=data.*s_r; % RC
clear s_r
% 产生RCMC信号
fr=d_fr*(0:Nr-1); % 距离向的频率采样点 (Hz)
s_rcmc=exp(-i*4*pi*f_dc*WaveLength/2*( (AzimuthTime)'*ones(1,Nr) )/c.*( ones(Na,1)*fr ));
clear fr
data=data.*s_rcmc; % RCMC
clear s_rcmc

% IFFT转为距离向上的时域信号
data=ifft(data,[],2);
%% (2) AC
data=fft(data,[],1); % 方位向上的FFT

% 产生方位向上的相关信号
% 方位向频斜率与距离相关性比较强时，uncomment以下三行
R=c*RangeTime/2+c*first_rg_cell/Fr/2+R0; % 快时间对应的斜距
K=-2*Vr^2/WaveLength./R; % 斜距对应的调频斜率
clear R
clear RangeTime AzimuthTime

% 方位向频斜率与距离相关性比较弱时，用固定的Ka进行AC
K=-Ka*ones(1,Nr);

% 计算方位向上匹配滤波器对应的时间
AzimuthMatchFilterTime=(0:Na_mf-1)/Fa; 
% 方位向上的LMF信号
s_a=exp(-i*2*pi*f_dc*((AzimuthMatchFilterTime-Ta/2)'*ones(1,Nr))+...
    i*pi*((AzimuthMatchFilterTime-Ta/2)'.^2)*K).*...
    (abs((AzimuthMatchFilterTime-Ta/2)'*ones(1,Nr))<=Ta/2); 
clear K AzimuthMatchFilterTime
% 补零至Na点长
s_a=[s_a;zeros(Na-Na_mf,Nr)]; 
% 方位向上LFM信号的频谱
s_a=fft(s_a,[],1); 
% 方位向上匹配滤波器的频谱
s_a=conj(s_a); 

data=data.*s_a; % AC
clear s_a;

% IFFT转为方位向上的时域信号
data=ifft(data,[],1);


%%
% ======================== The End ========================