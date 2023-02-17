clear;clc;close all;

M = 0.1:0.1:1.2;
PSD_gap = [10 11 13 9 10 14 13 12 12 10 8 9];
x = 0.1:0.01:1.2;
p1 =      -466.9;
p2 =        1757;
p3 =       -2128;
p4 =       476.8;
p5 =       849.3;
p6 =      -625.3;
p7 =       146.8;
p8 =      0.5455;
fitcurve = p1*x.^7 + p2*x.^6 + p3*x.^5 + p4*x.^4 + p5*x.^3 + ...
    p6*x.^2 + p7*x + p8;

figure(1)
plot(M,PSD_gap,'x')
xlabel('Modulation index');ylabel('PSD Gap [dB/Hz]');
axis([0.1 1.3 5 20])
set(gcf,'unit','centimeters','position',[10 10 17 11])
set(gca,'Position',[.12 .15 .8 .75]);
set(gca,'FontSize',12,'FontName','Arial');
hold on
plot(x,fitcurve)
legend('Real value','Fitting Curve')
