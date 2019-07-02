clc;
close all;
N=32;
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
%基于FFT的自相关函数快速计算方法
Uk=fft(un,2*N);      %对un进行2N点的FFT
Sk=(1/N)*abs(Uk).^2; %计算功率谱估计Sk
r0=ifft(Sk);         %对功率谱估计Sk求逆FFT
r1=[r0(N+2:2*N),r0(1:N)];%根据教材式（3.1.8）求得自相关函数
r=xcorr(un,N-1,'biased');%直接计算自相关函数

r1i=imag(r1);
r1r=real(r1); 
ri=imag(r);
rr=real(r);  
figure 
stem(-N+1:N-1,r1i);
title('基于FFT的自相关函数快速计算')
xlabel('m') 
ylabel('虚部') 
figure 
stem(-N+1:N-1,r1r); 
title('基于FFT的自相关函数快速计算') 
xlabel('m') 
ylabel('实部') 
figure 
stem(-N+1:N-1,ri); 
title('直接计算自相关函数?')
xlabel('m')
ylabel('虚部')
figure 
stem(-N+1:N-1,rr);
title('直接计算自相关函数?')
xlabel('m')
ylabel('实部')







