clc;
close all;
N=256;
noise=(randn(1,N)+randn(1,N))/sqrt(2);%产生零均值、方差为1的复高斯白噪声序列
f1=0.15;
f2=0.17;%信号的归一化频率
f3=0.26; 
SNR1=30;
SNR2=30; %信号的信噪比
SNR3=27;
A1=10^(SNR1/20);
A2=10^(SNR2/20); %信号的幅度
A3=10^(SNR3/20);
signal1=A1*exp(1i*2*pi*f1*(0:N-1));
signal2=A2*exp(1i*2*pi*f2*(0:N-1));%产生复正弦信号
signal3=A3*exp(1i*2*pi*f3*(0:N-1));
un=signal1+signal2+signal3+noise; %产生观察样本

%周期图法
NF=1024;%周期图法中FFT的点数
Spr=fftshift((1/NF)*abs(fft(un,NF)).^2);
Spr = abs(Spr);
Spr_dB = 10*log10(Spr/max(Spr));
w2 = linspace(-pi,pi,length(Spr_dB));
x2 = w2/2/pi;
figure;
plot(x2,Spr_dB);
xlim([-0.5 0.5])
xlabel('归一化频率/f');
ylabel('归一化功率谱/dB');
title('周期图法');

%BT法：自相关函数的单边长度M=64
M=64;        %自相关函数的单边长度
r=xcorr(un,M,'biased');%计算自相关函数
NF=1024;     %BT法中FFT的点数
BT=fftshift(fft(r,NF));%BT法计算功率谱
BT=abs(BT);
BT_dB = 10*log10(BT/max(BT));
w1 = linspace(-pi,pi,length(BT_dB));
x1 = w1/2/pi;
figure;
plot(x1,BT_dB);
xlim([-0.5 0.5]);
xlabel('归一化频率/f');
ylabel('归一化功率谱/dB');
title('BT法');