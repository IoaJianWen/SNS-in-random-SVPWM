n_u = length(u_A);
n_ufft = 32768;
u_fft = fft(u_A,n_ufft);
u_fft = abs(u_fft/n_u);
% y_fft = y_fft(1:(n+1)/2);
u_fft(2:end) = 2 * u_fft(2:end);
P_u = 20*log10(u_fft);
P_u_smooth = smoothdata(P_u,'gaussian',50);

n_i = length(i_A);
n_ifft = 32768;
i_fft = fft(i_A,n_ifft);
i_fft = abs(i_fft/n_i);
% y_fft = y_fft(1:(n+1)/2);
i_fft(2:end) = 2 * i_fft(2:end);
P_i = 20*log10(i_fft);
P_i_smooth = smoothdata(P_i,'gaussian',50);
% Y = 200 * Y;
% n_Y = length(Y);
% n_Yfft = 2048;
% RY = xcorr(Y);
% RY = 1/n_Y * RY(n_Y:end);
% p = fft(RY,n_Yfft);
% p = abs(p / n_Yfft);
% p(2:end) = 2 * p(2:end);
% P = 10 * log10(p);

T = 0:4/f0/n_u:(4/f0 - 1/f0/n_u);
f = 0:fss/n_ufft:(fss - fss/n_ufft);

% P = smoothdata(P_1,'gaussian',10);
%%
figure(1)
% subplot(2,1,1)
plot(T,u_A)
xlabel('Time [s]');ylabel('Amplitude [V]');
set(gcf,'unit','centimeters','position',[10 10 17 11])
set(gca,'Position',[.12 .15 .8 .75]);
set(gca,'FontSize',12,'FontName','Arial');

figure(2)
% subplot(2,1,2)
plot(f(1:(n_ufft/2+1)),P_u_smooth(1:(n_ufft/2+1)))
xlabel('Frequency [Hz]');ylabel('PSD [dB/Hz]');
axis([0 2e4 -60 15])
set(gcf,'unit','centimeters','position',[10 10 17 11])
set(gca,'Position',[.12 .15 .8 .75]);
set(gca,'FontSize',12,'FontName','Arial');

%%
figure(3)
% subplot(2,1,1)
plot(T,i_A)
hold on
plot(T,i_B)
hold on
plot(T,i_C)
xlabel('Time [s]');ylabel('Amplitude [A]');legend('i_a','i_b','i_c')
axis([0 0.08 -2.5 2.5])
set(gcf,'unit','centimeters','position',[10 10 17 11])
set(gca,'Position',[.12 .15 .8 .75]);
set(gca,'FontSize',12,'FontName','Arial');

figure(4)
% subplot(2,1,2)
plot(f(1:(n_ufft/2+1)),P_i_smooth(1:(n_ufft/2+1)))
xlabel('Frequency [Hz]');ylabel('PSD [dB/Hz]')
axis([0 2e4 -70 0])
set(gcf,'unit','centimeters','position',[10 10 17 11])
set(gca,'Position',[.12 .15 .8 .75]);
set(gca,'FontSize',12,'FontName','Arial');
% figure的position中的[left bottom width height] 是指figure的可画图的部分的左下角的坐标以及宽度和高度。