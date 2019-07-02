clc;
close all;
N=32;
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
%����FFT������غ������ټ��㷽��
Uk=fft(un,2*N);      %��un����2N���FFT
Sk=(1/N)*abs(Uk).^2; %���㹦���׹���Sk
r0=ifft(Sk);         %�Թ����׹���Sk����FFT
r1=[r0(N+2:2*N),r0(1:N)];%���ݽ̲�ʽ��3.1.8���������غ���
r=xcorr(un,N-1,'biased');%ֱ�Ӽ�������غ���

r1i=imag(r1);
r1r=real(r1); 
ri=imag(r);
rr=real(r);  
figure 
stem(-N+1:N-1,r1i);
title('����FFT������غ������ټ���')
xlabel('m') 
ylabel('�鲿') 
figure 
stem(-N+1:N-1,r1r); 
title('����FFT������غ������ټ���') 
xlabel('m') 
ylabel('ʵ��') 
figure 
stem(-N+1:N-1,ri); 
title('ֱ�Ӽ�������غ���?')
xlabel('m')
ylabel('�鲿')
figure 
stem(-N+1:N-1,rr);
title('ֱ�Ӽ�������غ���?')
xlabel('m')
ylabel('ʵ��')







