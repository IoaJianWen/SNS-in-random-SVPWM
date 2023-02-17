n_u_A = length(u_A);
n_u_gap = length(u_A_gap);

n_ufft = 32768;

u_A_fft = fft(u_A,n_ufft);
u_A_fft = abs(u_A_fft/n_u_A);
% y_fft = y_fft(1:(n+1)/2);
u_A_fft(2:end) = 2 * u_A_fft(2:end);
P_u_A = 20*log10(u_A_fft);
P_u_A = P_u_A(1:(n_ufft/2 + 1));
P_u_A_smooth = smoothdata(P_u_A,'gaussian',50);

u_A_gap_fft = fft(u_A_gap,n_ufft);
u_A_gap_fft = abs(u_A_gap_fft/n_u_gap);
% y_fft = y_fft(1:(n+1)/2);
u_A_gap_fft(2:end) = 2 * u_A_gap_fft(2:end);
P_u_A_gap = 20*log10(u_A_gap_fft);
P_u_A_gap = P_u_A_gap(1:(n_ufft/2 + 1));
P_u_A_gap_smooth = smoothdata(P_u_A_gap,'gaussian',50);

n_i_A = length(i_A);

n_ifft = 32768;

i_A_fft = fft(i_A,n_ifft);
i_A_fft = abs(i_A_fft/n_i_A);
i_A_fft(2:end) = 2 * i_A_fft(2:end);
P_i_A = 20*log10(i_A_fft);
P_i_A = P_i_A(1:(n_ufft/2 + 1));
P_i_A_smooth = smoothdata(P_i_A,'gaussian',50);

n_i_A_gap = length(i_A_gap);

i_A_gap_fft = fft(i_A_gap,n_ifft);
i_A_gap_fft = abs(i_A_gap_fft/n_i_A_gap);
% y_fft = y_fft(1:(n+1)/2);
i_A_gap_fft(2:end) = 2 * i_A_gap_fft(2:end);
P_i_A_gap = 20*log10(i_A_gap_fft);
P_i_A_gap = P_i_A_gap(1:(n_ifft/2 + 1));
P_i_A_gap_smooth = smoothdata(P_i_A_gap,'gaussian',50);

% Y = 200 * Y;
% n_Y = length(Y);
% n_Yfft = 2048;
% RY = xcorr(Y);
% RY = 1/n_Y * RY(n_Y:end);
% p = fft(RY,n_Yfft);
% p = abs(p / n_Yfft);
% p(2:end) = 2 * p(2:end);
% P = 10 * log10(p);

T_A = 0:4/f0/n_u_A:(4/f0 - 1/f0/n_u_A);
f = 0:fss/n_ufft:(fss - fss/n_ufft);
f = f(1:(n_ufft/2 + 1));


%%
figure(1)
plot(f,P_u_A_smooth,'--','LineWidth',0.8)
xlabel('Frequency [Hz]');ylabel('PSD [dB/Hz]');
axis([0 12e3 -35 15])
set(gcf,'unit','centimeters','position',[10 10 17 11])
set(gca,'Position',[.12 .15 .8 .75]);
set(gca,'FontSize',12,'FontName','Arial');
hold on
plot(f,P_u_A_gap_smooth)
legend('RP-SVPWM','the first SNS method')

%%
figure(2)
plot(f,P_i_A_smooth)
xlabel('Frequency [Hz]');ylabel('PSD [dB/Hz]');
axis([0 12e3 -80 -10])
set(gcf,'unit','centimeters','position',[10 10 17 11])
set(gca,'Position',[.12 .15 .8 .75]);
set(gca,'FontSize',12,'FontName','Arial');
hold on
plot(f,P_i_A_gap_smooth)
legend('RP-SVPWM','the first SNS method')

%%
% figure(1)
% % subplot(2,1,1)
% plot(T,u_A)
% xlabel('Time [s]');ylabel('Amplitude [V]');
% set(gcf,'unit','centimeters','position',[10 10 17 11])
% set(gca,'Position',[.12 .15 .8 .75]);
% set(gca,'FontSize',12,'FontName','Arial');

figure(2)
% subplot(2,1,2)
plot(f(1:8193),P_u_A_gap(1:8193)) % 0-10k
xlabel('Frequency [Hz]');ylabel('PSD [dB/Hz]');
axis([0 10e3 -60 30])
set(gcf,'unit','centimeters','position',[10 10 17 11])
set(gca,'Position',[.12 .15 .8 .75]);
set(gca,'FontSize',12,'FontName','Arial');

% figure(3)
% plot(f(4915:6555),P_u(4915:6555)) % 6k-8k
% xlabel('Frequency [Hz]');ylabel('PSD [dB/Hz]');
% axis([6e3 8e3 -50 0])
% set(gcf,'unit','centimeters','position',[10 10 17 11])
% set(gca,'Position',[.12 .15 .8 .75]);
% set(gca,'FontSize',12,'FontName','Arial');
% 
% figure(4)
% plot(f(3277:8193),P_u(3277:8193)) % 4k-10k
% xlabel('Frequency [Hz]');ylabel('PSD [dB/Hz]');
% axis([4e3 10e3 -50 0])
% set(gcf,'unit','centimeters','position',[10 10 17 11])
% set(gca,'Position',[.12 .15 .8 .75]);
% set(gca,'FontSize',12,'FontName','Arial');
%%
figure(5)
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

figure(6)
% subplot(2,1,2)
plot(f(1:(n_ufft/2+1)),P_i(1:(n_ufft/2+1)))
xlabel('Frequency [Hz]');ylabel('PSD [dB/Hz]')
axis([0 2e4 -100 10])
set(gcf,'unit','centimeters','position',[10 10 17 11])
set(gca,'Position',[.12 .15 .8 .75]);
set(gca,'FontSize',12,'FontName','Arial');
% figure的position中的[left bottom width height] 是指figure的可画图的部分的左下角的坐标以及宽度和高度。