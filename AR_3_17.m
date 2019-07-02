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

p=16;                     %ARģ�͵Ľ���
k = zeros(1,p);           % ��ʼ������ϵ��k(m)
Sigma= zeros(1,p);        % ��ʼ��p��ARģ�͵���������

r0 = xcorr(un,p,'biased');%ֱ�Ӽ�������غ���
r = r0(p+1:2*p+1);        %��ȡr(0),r(1),������,r(p)
%����һ��ARģ�͵�ϵ�������뷽��
a(1,1) = -r(2)/r(1);      %1��ARģ�͵�ϵ��
Sigma(1) = r(1)-(abs(r(2))^2)/r(1);%1��ARģ�͵����뷽��

% Levinson-Durbin�����㷨
 for m=2:p
    k(m) = -(r(m+1) + sum(a(m-1,1:m-1).*r(m:-1:2)))/Sigma(m-1);
    a(m,m) = k(m);
    for i = 1 : m-1
        a(m,i) = a(m-1, i)+k(m)*conj(a(m-1,m-i));
     end
   Sigma(m) = Sigma(m-1)*(1-abs(k(m))^2);
 end

 %����16��ARģ�͵Ĺ�����
 NF = 1024;  %AR������FFT�ĵ���
 Par = Sigma(p)./abs(fftshift(fft([1 a(p,:)],NF))).^2;%P��ARģ�͵Ĺ�����
 Par = abs(Par);
 Par_dB = 10*log10(Par/max(Par));
 w = linspace(-pi,pi,NF);
 x = w/2/pi;
 figure;
 plot(x,Par_dB);
 xlim([-0.5 0.5]);
 xlabel('��һ��Ƶ��/f');
 ylabel('��һ��������/dB');
 title('16��ARģ�͵Ĺ����׹���');