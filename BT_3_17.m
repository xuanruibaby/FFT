clc;
close all;
N=256;
noise=(randn(1,N)+randn(1,N))/sqrt(2);%�������ֵ������Ϊ1�ĸ���˹����������
f1=0.15;
f2=0.17;%�źŵĹ�һ��Ƶ��
f3=0.26; 
SNR1=30;
SNR2=30; %�źŵ������
SNR3=27;
A1=10^(SNR1/20);
A2=10^(SNR2/20); %�źŵķ���
A3=10^(SNR3/20);
signal1=A1*exp(1i*2*pi*f1*(0:N-1));
signal2=A2*exp(1i*2*pi*f2*(0:N-1));%�����������ź�
signal3=A3*exp(1i*2*pi*f3*(0:N-1));
un=signal1+signal2+signal3+noise; %�����۲�����

%����ͼ��
NF=1024;%����ͼ����FFT�ĵ���
Spr=fftshift((1/NF)*abs(fft(un,NF)).^2);
Spr = abs(Spr);
Spr_dB = 10*log10(Spr/max(Spr));
w2 = linspace(-pi,pi,length(Spr_dB));
x2 = w2/2/pi;
figure;
plot(x2,Spr_dB);
xlim([-0.5 0.5])
xlabel('��һ��Ƶ��/f');
ylabel('��һ��������/dB');
title('����ͼ��');

%BT��������غ����ĵ��߳���M=64
M=64;        %����غ����ĵ��߳���
r=xcorr(un,M,'biased');%��������غ���
NF=1024;     %BT����FFT�ĵ���
BT=fftshift(fft(r,NF));%BT�����㹦����
BT=abs(BT);
BT_dB = 10*log10(BT/max(BT));
w1 = linspace(-pi,pi,length(BT_dB));
x1 = w1/2/pi;
figure;
plot(x1,BT_dB);
xlim([-0.5 0.5]);
xlabel('��һ��Ƶ��/f');
ylabel('��һ��������/dB');
title('BT��');