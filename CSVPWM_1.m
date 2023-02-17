%% 
% 传统SVPWM 线电压信号波形与功率谱 三相电流信号波形与功率谱

clear;clc;close all;

fss = 40e3; % 离散采样率
f0 = 50; % 合成旋转电压的频率
fs = 2500; % 固定开关频率
M = 0.7;

% A相开关信号
T_duty_A = 0;
duty_A = [];
S_A = [];
ite_A = 1;
while (T_duty_A < 4/f0)
    if (T_duty_A >= 0 && T_duty_A < 1/(6*f0)) || ...
            (T_duty_A >= 6/(6*f0) && T_duty_A < 7/(6*f0)) ||...
            (T_duty_A >= 12/(6*f0) && T_duty_A < 13/(6*f0)) ||...
            (T_duty_A >= 18/(6*f0) && T_duty_A < 19/(6*f0))
        duty_A(ite_A) = 1/2 * (1 + M * sin(2 * pi * f0 * T_duty_A + pi / 3));

    elseif (T_duty_A >= 1/(6*f0) && T_duty_A < 2/(6*f0)) ||...
            (T_duty_A >= 7/(6*f0) && T_duty_A < 8/(6*f0)) ||...
            (T_duty_A >= 13/(6*f0) && T_duty_A < 14/(6*f0)) ||...
            (T_duty_A >= 19/(6*f0) && T_duty_A < 20/(6*f0))
        duty_A(ite_A) = 1/2 * (1 + M * sqrt(3) * cos(2 * pi * f0 * T_duty_A));

    elseif (T_duty_A >= 2/(6*f0) && T_duty_A < 3/(6*f0)) ||...
            (T_duty_A >= 8/(6*f0) && T_duty_A < 9/(6*f0)) ||...
            (T_duty_A >= 14/(6*f0) && T_duty_A < 15/(6*f0)) ||...
            (T_duty_A >= 20/(6*f0) && T_duty_A < 21/(6*f0))
        duty_A(ite_A) = 1/2 * (1 + M * sin(pi / 3 - 2 * pi * f0 * T_duty_A));

    elseif (T_duty_A >= 3/(6*f0) && T_duty_A < 4/(6*f0)) ||...
            (T_duty_A >= 9/(6*f0) && T_duty_A < 10/(6*f0)) ||...
            (T_duty_A >= 15/(6*f0) && T_duty_A < 16/(6*f0)) ||...
            (T_duty_A >= 21/(6*f0) && T_duty_A < 22/(6*f0))
        duty_A(ite_A) = 1/2 * (1 + M * sin(2 * pi * f0 * T_duty_A + pi / 3));

    elseif (T_duty_A >= 4/(6*f0) && T_duty_A < 5/(6*f0)) ||...
            (T_duty_A >= 10/(6*f0) && T_duty_A < 11/(6*f0)) ||...
            (T_duty_A >= 16/(6*f0) && T_duty_A < 17/(6*f0)) ||...
            (T_duty_A >= 22/(6*f0) && T_duty_A < 23/(6*f0))
        duty_A(ite_A) = 1/2 * (1 + M * sqrt(3) * cos(2 * pi * f0 * T_duty_A));

    else
        duty_A(ite_A) = 1/2 * (1 + M * sin(pi / 3 - 2 * pi * f0 * T_duty_A));
    end
    
    pulsewidth = (1/fs) * duty_A(ite_A);
    
    T_ss = 0:1/fss:(1/fs - 1/fss);
    rect = (0:1/fs:1/fs)'; % 一次循环生成一个 该时刻占空比和fs 对应下的脉冲
    Y_1 = pulstran(T_ss - 1/(2*fs),rect,'rectpuls',pulsewidth);
    S_A = [S_A,Y_1];
    
    T_duty_A = T_duty_A + 1/fs; % 占空比的采样间隔（时间）
    ite_A = ite_A + 1;
end

% B相开关信号
T_duty_B = 0;
duty_B = [];
S_B = [];
ite_B = 1;
while (T_duty_B < 4/f0)
    if (T_duty_B >= 0 && T_duty_B < 1/(6*f0)) || ...
            (T_duty_B >= 6/(6*f0) && T_duty_B < 7/(6*f0)) ||...
            (T_duty_B >= 12/(6*f0) && T_duty_B < 13/(6*f0)) ||...
            (T_duty_B >= 18/(6*f0) && T_duty_B < 19/(6*f0))
        duty_B(ite_B) = 1/2 * (1 + M * sqrt(3) * sin(2 * pi * f0 * T_duty_B - pi / 6));

    elseif (T_duty_B >= 1/(6*f0) && T_duty_B < 2/(6*f0)) ||...
            (T_duty_B >= 7/(6*f0) && T_duty_B < 8/(6*f0)) ||...
            (T_duty_B >= 13/(6*f0) && T_duty_B < 14/(6*f0)) ||...
            (T_duty_B >= 19/(6*f0) && T_duty_B < 20/(6*f0))
        duty_B(ite_B) = 1/2 * (1 + M * sin(2 * pi * f0 * T_duty_B));

    elseif (T_duty_B >= 2/(6*f0) && T_duty_B < 3/(6*f0)) ||...
            (T_duty_B >= 8/(6*f0) && T_duty_B < 9/(6*f0)) ||...
            (T_duty_B >= 14/(6*f0) && T_duty_B < 15/(6*f0)) ||...
            (T_duty_B >= 20/(6*f0) && T_duty_B < 21/(6*f0))
        duty_B(ite_B) = 1/2 * (1 + M * sin(2 * pi * f0 * T_duty_B - pi / 3));

    elseif (T_duty_B >= 3/(6*f0) && T_duty_B < 4/(6*f0)) ||...
            (T_duty_B >= 9/(6*f0) && T_duty_B < 10/(6*f0)) ||...
            (T_duty_B >= 15/(6*f0) && T_duty_B < 16/(6*f0)) ||...
            (T_duty_B >= 21/(6*f0) && T_duty_B < 22/(6*f0))
        duty_B(ite_B) = 1/2 * (1 + M * sqrt(3) * sin(2 * pi * f0 * T_duty_B - pi / 6));

    elseif (T_duty_B >= 4/(6*f0) && T_duty_B < 5/(6*f0)) ||...
            (T_duty_B >= 10/(6*f0) && T_duty_B < 11/(6*f0)) ||...
            (T_duty_B >= 16/(6*f0) && T_duty_B < 17/(6*f0)) ||...
            (T_duty_B >= 22/(6*f0) && T_duty_B < 23/(6*f0))
        duty_B(ite_B) = 1/2 * (1 + M * sin(2 * pi * f0 * T_duty_B));

    else
        duty_B(ite_B) = 1/2 * (1 + M * sin(2 * pi * f0 * T_duty_B - pi / 3));
    end
    
    pulsewidth = (1/fs) * duty_B(ite_B);
    
    T_ss = 0:1/fss:(1/fs - 1/fss);
    rect = (0:1/fs:1/fs)'; % 一次循环生成一个 该时刻占空比和fs 对应下的脉冲
    Y_2 = pulstran(T_ss - 1/(2*fs),rect,'rectpuls',pulsewidth);
    S_B = [S_B,Y_2];
    
    T_duty_B = T_duty_B + 1/fs; % 占空比的采样间隔（时间）
    ite_B = ite_B + 1;
end

% C相开关信号
T_duty_C = 0;
duty_C = [];
S_C = [];
ite_C = 1;
while (T_duty_C < 4/f0)
    if (T_duty_C >= 0 && T_duty_C < 1/(6*f0)) || ...
            (T_duty_C >= 6/(6*f0) && T_duty_C < 7/(6*f0)) ||...
            (T_duty_C >= 12/(6*f0) && T_duty_C < 13/(6*f0)) ||...
            (T_duty_C >= 18/(6*f0) && T_duty_C < 19/(6*f0))
        duty_C(ite_C) = 1/2 * (1 - M * sin(2 * pi * f0 * T_duty_C + pi / 3));

    elseif (T_duty_C >= 1/(6*f0) && T_duty_C < 2/(6*f0)) ||...
            (T_duty_C >= 7/(6*f0) && T_duty_C < 8/(6*f0)) ||...
            (T_duty_C >= 13/(6*f0) && T_duty_C < 14/(6*f0)) ||...
            (T_duty_C >= 19/(6*f0) && T_duty_C < 20/(6*f0))
        duty_C(ite_C) = 1/2 * (1 - M * sin(2 * pi * f0 * T_duty_C));

    elseif (T_duty_C >= 2/(6*f0) && T_duty_C < 3/(6*f0)) ||...
            (T_duty_C >= 8/(6*f0) && T_duty_C < 9/(6*f0)) ||...
            (T_duty_C >= 14/(6*f0) && T_duty_C < 15/(6*f0)) ||...
            (T_duty_C >= 20/(6*f0) && T_duty_C < 21/(6*f0))
        duty_C(ite_C) = 1/2 * (1 - M * sqrt(3) * sin(2 * pi * f0 * T_duty_C + pi / 6));

    elseif (T_duty_C >= 3/(6*f0) && T_duty_C < 4/(6*f0)) ||...
            (T_duty_C >= 9/(6*f0) && T_duty_C < 10/(6*f0)) ||...
            (T_duty_C >= 15/(6*f0) && T_duty_C < 16/(6*f0)) ||...
            (T_duty_C >= 21/(6*f0) && T_duty_C < 22/(6*f0))
        duty_C(ite_C) = 1/2 * (1 - M * sin(2 * pi * f0 * T_duty_C + pi / 3));

    elseif (T_duty_C >= 4/(6*f0) && T_duty_C < 5/(6*f0)) ||...
            (T_duty_C >= 10/(6*f0) && T_duty_C < 11/(6*f0)) ||...
            (T_duty_C >= 16/(6*f0) && T_duty_C < 17/(6*f0)) ||...
            (T_duty_C >= 22/(6*f0) && T_duty_C < 23/(6*f0))
        duty_C(ite_C) = 1/2 * (1 - M * sin(2 * pi * f0 * T_duty_C));

    else
        duty_C(ite_C) = 1/2 * (1 - M * sqrt(3) * sin(2 * pi * f0 * T_duty_C + pi / 6));
    end
    
    pulsewidth = (1/fs) * duty_C(ite_C);
    
    T_ss = 0:1/fss:(1/fs - 1/fss);
    rect = (0:1/fs:1/fs)'; % 一次循环生成一个 该时刻占空比和fs 对应下的脉冲
    Y_3 = pulstran(T_ss - 1/(2*fs),rect,'rectpuls',pulsewidth);
    S_C = [S_C,Y_3];
    
    T_duty_C = T_duty_C + 1/fs; % 占空比的采样间隔（时间）
    ite_C = ite_C + 1;
end

% 相电压 u_A u_B u_C
u_dc = 24;
u_A = u_dc * (1/3) * (2 * S_A - S_B - S_C);
u_B = u_dc * (1/3) * (2 * S_B - S_A - S_C);
u_C = u_dc * (1/3) * (2 * S_C - S_A - S_B);

% 相电流 i_a i_b i_c
L = 5e-3;
R = 5;
b = 1/L/fss;
a = [1,-exp(-R/L*(1/fss))];
i_A = filter(b,a,u_A);
i_B = filter(b,a,u_B);
i_C = filter(b,a,u_C);
