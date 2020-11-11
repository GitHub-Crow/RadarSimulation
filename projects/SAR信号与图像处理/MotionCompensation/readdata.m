%% ��ȡ����
clc; close all;
C=3e8;                   % ����
Vr=154.195864;           % �����ٶ�
Tr=2.4e-6;               % ������
lambda=0.03125;          % �ز�����
f0=C/lambda;             % �ز�Ƶ�ʣ�X����
Kr=-2e14;                % �������Ƶ��
Bw=Kr*Tr;                % ������LFM����
Fr=548571428.571429;     % �����������
PRF=533.330793;          % ��λ�������
theta=0.04;              % ��λ�����Ƕȣ�����
sin_theta=sin(theta/2);  
range_size=16384;
azimuth_size=20480;
Nr=range_size;           % �������������
Na=azimuth_size/8;       % ��λ���������
near_range=23306.25;     % ��һ�����������
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
forward = m(4:12:end);          % ǰ������
cross = m(5:12:end);            % ����
height = m(6:12:end);           % �߶���
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
