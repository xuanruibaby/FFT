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

p=16;                     %AR模型的阶数
k = zeros(1,p);           % 初始化反射系数k(m)
Sigma= zeros(1,p);        % 初始化p阶AR模型的噪声方差

r0 = xcorr(un,p,'biased');%直接计算自相关函数
r = r0(p+1:2*p+1);        %提取r(0),r(1),···,r(p)
%计算一阶AR模型的系数与输入方差
a(1,1) = -r(2)/r(1);      %1阶AR模型的系数
Sigma(1) = r(1)-(abs(r(2))^2)/r(1);%1阶AR模型的输入方差

% Levinson-Durbin迭代算法
 for m=2:p
    k(m) = -(r(m+1) + sum(a(m-1,1:m-1).*r(m:-1:2)))/Sigma(m-1);
    a(m,m) = k(m);
    for i = 1 : m-1
        a(m,i) = a(m-1, i)+k(m)*conj(a(m-1,m-i));
     end
   Sigma(m) = Sigma(m-1)*(1-abs(k(m))^2);
 end

 %计算16阶AR模型的功率谱
 NF = 1024;  %AR方法中FFT的点数
 Par = Sigma(p)./abs(fftshift(fft([1 a(p,:)],NF))).^2;%P阶AR模型的功率谱
 Par = abs(Par);
 Par_dB = 10*log10(Par/max(Par));
 w = linspace(-pi,pi,NF);
 x = w/2/pi;
 figure;
 plot(x,Par_dB);
 xlim([-0.5 0.5]);
 xlabel('归一化频率/f');
 ylabel('归一化功率谱/dB');
 title('16阶AR模型的功率谱估计');