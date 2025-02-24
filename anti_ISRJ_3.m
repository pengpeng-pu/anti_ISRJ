%————————————波形设计————————————————————————————
% 清除工作空间和命令窗口
clear;
clc;
%% 信号参数
T = 40e-6; % 信号时宽 40 us
B = 40e6; % 带宽 40 MHz
fs = 2 * B; % 采样频率，至少为带宽的两倍
t = 0 : 1/fs : T-1/fs; % 时间向量
L = length(t); % 信号长度（采样点）
% 距离参数
R_target1 = 2e3; % 目标1距离中心 2 km
R_target2 = 1.5e3;
R_jammer = 1.5e3; % 干扰机距离中心 1.5 km
c = 3e8; % 光速
target1_delay = 2 * R_target1 / c; % 往返延迟
target2_delay = 2 * R_target2 / c; % 往返延迟
jammer_delay = 2 * R_jammer / c;
jammer_delay_index = round(jammer_delay * fs);%干扰机路程时延
% 干扰机参数
T_samp_jammer = 2e-6; % 干扰机采样时宽 2 us
Ts_jammer = 8e-6; % 干扰机采样周期 8 us
% 信干比和信噪比
SIR = -15; % 信干比 -15 dB
SNR = 0; % 信噪比 0 dB
RCS = 1;
% 生成间歇采样脉冲串与干扰信号
pulse_train = zeros(length(t),1);
t_samp_jammer = 0 : 1/fs : T_samp_jammer-1/fs;
num_samps_jammer = length(t_samp_jammer);% 单次采样点数
num_periods = floor((T - T_samp_jammer) / Ts_jammer)+1;% 信号时宽内采样次数
for n = 0:num_periods-1%第n次采样
    start_idx = round((n * Ts_jammer) * fs)+1;%单次采样起始采样点数
    end_idx = start_idx + num_samps_jammer - 1;%终止采样点数
    if end_idx > L
        end_idx = L;
    end
    pulse_train(start_idx:end_idx) = 1; % 将采样时段置为 1
end
f0 = 10e6; % 起始频率
k = B / T; % 调频斜率
LFM_signal = exp(1j * (2 * pi * f0 * t + pi * k * t.^2));% LFM 信号
LFM_signal = LFM_signal';

signal_power = sum(abs(LFM_signal).^2)/L;% 计算信号功率
interference_power = signal_power / (10^(SIR/20));% 计算干扰功率
noise_power = signal_power / (10^(SNR/20));% 计算噪声功率
LFM_j = LFM_signal .* pulse_train;% 3200*1
Jammer_signal = LFM_j * sqrt(interference_power/sum(abs(LFM_j).^2)*L);
%% ISCRJ
jammer_ISCRJ = zeros(L,1);% 3200*1
Jammer_ISCRJ = zeros(L,1);% 3200*1
Jammer_n = zeros(5*num_samps_jammer,1);
z = length(Jammer_n);
for n = 0:num_periods-1% 第n次采样
    start_idx = round((n * Ts_jammer) * fs)+1;
    end_idx = start_idx + num_samps_jammer-1;
    Jammer_n = [Jammer_signal(start_idx : end_idx,1);Jammer_n(1:z-end_idx+start_idx-1,1)];
    if end_idx+1+z-1 <= L
        Jammer_ISCRJ(end_idx+1 : end_idx+1+z-1) = Jammer_ISCRJ(end_idx+1 : end_idx+1+z-1) + Jammer_n;
    else
        Jammer_ISCRJ(end_idx+1 : L,1) = Jammer_ISCRJ(end_idx+1 : L,1) + Jammer_n(1:L-end_idx,1);
    end
end
if jammer_delay_index < length(t)
    jammer_ISCRJ(jammer_delay_index:end,1) = Jammer_ISCRJ(1:end-jammer_delay_index + 1,1);
end
% 目标回波(多目标)
target1_echo = zeros(L,1);% 3200*1
target2_echo = zeros(L,1);% 3200*1
target1_delay_index = round(target1_delay * fs);
target2_delay_index = round(target2_delay * fs);
if target1_delay_index < length(t) && target2_delay_index < length(t)
    target1_echo(target1_delay_index:end,1) = RCS * LFM_signal(1:end-target1_delay_index + 1,1);%目标1回波
    target2_echo(target2_delay_index:end,1) = RCS * LFM_signal(1:end-target2_delay_index + 1,1);%目标2回波
end
% 计算噪声信号
noise = sqrt(noise_power/2)*(randn(L,1) + 1j * randn(L,1));
sig1 = target1_echo+jammer_ISCRJ+noise;
sig2 = target1_echo+target2_echo+jammer_ISCRJ+noise;
% 匹配滤波
matching_filter = conj(flip(LFM_signal));
output_single_ISCRJ = conv(sig1, matching_filter);
output_single_ISCRJ = fftshift(output_single_ISCRJ);
figure;
x = 1:length(output_single_ISCRJ);
plot(2*x,abs(output_single_ISCRJ)/max(abs(output_single_ISCRJ)));
xlim([0,6000]);
title('LFM信号匹配滤波输出');
xlabel('距离(m)');
ylabel('幅度(dB)');

output_single_ISCRJ = conv(sig2, matching_filter);
output_single_ISCRJ = fftshift(output_single_ISCRJ);
figure;
plot(2*x,abs(output_single_ISCRJ)/max(abs(output_single_ISCRJ)));
xlim([0,6000]);
title('LFM信号匹配滤波输出');
xlabel('距离(m)');
ylabel('幅度(dB)');

% 非匹配滤波
non_matched_filter_impulse = conj(flip(LFM_signal - LFM_j));
filtered_multi_ISCRJ = conv(sig1, non_matched_filter_impulse);
filtered_multi_ISCRJ = fftshift(filtered_multi_ISCRJ);
figure;
plot(2*x,abs(filtered_multi_ISCRJ)/max(abs(filtered_multi_ISCRJ)));
xlim([0,6000]);
title('LFM信号非匹配滤波输出');
xlabel('距离(m)');
ylabel('幅度(dB)');

filtered_multi_ISCRJ = conv(sig2, non_matched_filter_impulse);
filtered_multi_ISCRJ = fftshift(filtered_multi_ISCRJ);
figure;
plot(2*x,abs(filtered_multi_ISCRJ)/max(abs(filtered_multi_ISCRJ)));
xlim([0,6000]);
title('LFM信号非匹配滤波输出');
xlabel('距离(m)');
ylabel('幅度(dB)');