%% 读取数据
clc; close all;
C=3e8;                   % 光速
Vr=154.195864;           % 飞行速度
Tr=2.4e-6;               % 脉冲宽度
lambda=0.03125;          % 载波波长
f0=C/lambda;             % 载波频率，X波段
Kr=-2e14;                % 距离向调频率
Bw=Kr*Tr;                % 距离向LFM带宽
Fr=548571428.571429;     % 距离向采样率
PRF=533.330793;          % 方位向采样率
theta=0.04;              % 方位波束角度，弧度
sin_theta=sin(theta/2);  
range_size=16384;
azimuth_size=20480;
Nr=range_size;           % 距离向采样点数
Na=azimuth_size/8;       % 方位向采样点数
near_range=23306.25;     % 第一个采样点距离
hang=0;
lie=1;

fid=fopen('data_before_moco.dat','rb');
data=zeros(Na,Nr);
%data=zeros(Na,Nr+2*nr);
fseek(fid,hang*Na*Nr*8,-1);
for index=1:Na
    m=fread(fid,[1,Nr*2],'float32');
    data(index,:)=m(1:2:end)+1j*m(2:2:end);
end
fclose(fid);
figure;imagesc(abs(data));  
fid = fopen('mocodata.dat','rb');
head = fread(fid,[1,9],'double');
m = fread(fid,[1,12*20480],'double');
time = m(1:12:end);
forward = m(4:12:end);          % 前进方向
cross = m(5:12:end);            % 横向
height = m(6:12:end);           % 高度向
ref_cross = zeros(1,20480);
ref_height = m(10:12:end);
ref_forward = m(7:12:end);
time = time(1 : Na);
forward = forward(1 : Na);
cross = cross(1 : Na);
height = height(1 : Na);
ref_cross = ref_cross(1 : Na);
ref_height = ref_height(1 : Na);
ref_forward = ref_forward(1 : Na);
plot3(forward, cross, height);
hold on;
plot3(ref_forward, ref_cross, ref_height,'r');
xlabel('xx');
ylabel('yy');
grid on;
fclose(fid);
