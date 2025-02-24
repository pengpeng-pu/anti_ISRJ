%————————————波形设计————————————————————————————
% 清除工作空间和命令窗口
clear;
clc;
%% 信号参数
T = 40e-6; % 信号时宽 40 us
B = 40e6; % 带宽 40 MHz
fs = 2 * B; % 采样频率，至少为带宽的两倍
t = 0 : 1/fs : T-1/fs; % 时间向量
t1 = -T+1/fs : 1/fs : T-1/fs;
L = length(t); % 信号长度（采样点）
RCS = 1;
% 距离参数
R_center = 60e3; % 场景中心距离 60 km
R_target1 = 2e3; % 目标1距离中心 2 km
R_target2 = 1.5e3;% 目标2距离中心 1.5 km
R_jammer = 1.5e3; % 干扰机距离中心 1.5 km
c = 3e8; % 光速
target1_delay = 2 * R_target1 / c; % 往返延迟
target2_delay = 2 * R_target2 / c; % 往返延迟
jammer_delay = 2 * R_jammer / c;
% 干扰机参数
T_samp_jammer = 2e-6; % 干扰机采样时宽 2 us
Ts_jammer = 8e-6; % 干扰机采样周期 8 us
% 信干比和信噪比
SIR = -15; % 信干比 -15 dB
SNR = 0; % 信噪比 0 dB
%% 生成间歇采样脉冲串与干扰信号
pulse_train = zeros(size(t));
t_samp_jammer = 0 : 1/fs : T_samp_jammer-1/fs;
num_samps_jammer = length(t_samp_jammer);
num_periods = floor((T - T_samp_jammer) / Ts_jammer)+1;% 周期内采样次数
for n = 0:num_periods-1 % 第n次采样
    start_idx = round((n * Ts_jammer) * fs)+1;% 单次采样起始采样点数
    end_idx = start_idx + num_samps_jammer - 1;% 终止采样点数
    if end_idx > L
        end_idx = L;
    end
    pulse_train(start_idx:end_idx) = 1; % 将采样时段置为 1
end
pulse_train = pulse_train';
%% 迭代算法实现联合设计
% 利用随机相位编码信号初始化x,y,h
rng(42);
x = exp(1j * 2 * pi * rand(L, 1)); % 初始发射波形
h = exp(1j * 2 * pi * rand(L, 1)); % 初始非匹配滤波器
y = x;
iter_max = 100; % 最大迭代次数
b_max = 1;         % 脉冲压缩归一化峰值 0 dB对应线性值1
b_min = 10^(-30/20); % 干扰归一化峰值 -30 dB对应约0.0316
delta = 1e-6; % 截止条件
prev_total = 0;
Integrate_sidelobe = zeros(1,iter_max+1);
Integrate_jamming = zeros(1,iter_max+1);
H_SL = Matrix_HSL(L, h);
H_j = Matrix_HJ(L, pulse_train, h);
% 目标信号积分旁瓣
r = H_SL*x;
Integrate_sidelobe(1,1) = r' * r;
% 干扰信号积分能量
r_j = H_j*x;
Integrate_jamming(1,1) = r_j' * r_j;
Converg_Iter = 0;
steps = iter_max;
a = waitbar(0,'计算进度');

for iter = 1:iter_max
    iter
    rho = 1.2 ^ (iter-1); % 步长
    if rho > 500
        rho = 500;
    end
    % 构建矩阵和函数【步骤2】
    x_j = x .* pulse_train;
    X_SL = Matrix_XSL(L, x);
    X_j = Matrix_XJ(L, pulse_train, x);
    Q = X_SL' * X_SL + X_j' * X_j;
    % 提前计算 Q 的伪逆
    Q_inv = pinv(Q + 1e-10 * eye(size(Q)));
    % 固定 x，更新 h
    h = mismatched_filter_design(x , x_j , Q_inv , b_max , b_min);
    H_SL = Matrix_HSL(L, h);
    H_j = Matrix_HJ(L, pulse_train, h);
    h_p = h;
    h_j_p = h_p .* pulse_train;
    % 提前计算部分矩阵乘法
    H_SLtH_SL = H_SL' * H_SL;
    H_jtH_j = H_j' * H_j;
    h_p_h_pT = h_p * h_p';
    h_j_p_h_j_pT = h_j_p * h_j_p'; 
    % 计算分子 1
    numerator1 = 2 * H_SLtH_SL + 2 * H_jtH_j + rho * eye(L) + rho * h_p_h_pT + rho * h_j_p_h_j_pT;
    % 为避免数值不稳定，添加一个小的正则化项
    numerator1 = numerator1 + 1e-10 * eye(size(numerator1));
    x_numerator1 = pinv(numerator1);
    % 初始化拉格朗日对偶变量
    lambda = zeros(L, 1); 
    mu = 0;
    nu = 0;
    condition = false;
    while ~condition
        % 固定 h，更新 x【步骤3】
        x_prev = x;
        x = waveform_design(x_numerator1, h_p, h_j_p, x_prev, lambda, mu, nu, rho, b_max, b_min);
        % 更新中间变量 y【步骤4】
        z = -(x + lambda);
        y = -z ./ abs(z);
        % 更新拉格朗日对偶变量【步骤5】
        lambda = lambda + x - y;
        mu = mu + h_p' * x - b_max;
        nu = nu + h_j_p' * x - b_min;
        % 截止条件【步骤6】
        if norm(x - y)^2 <= delta 
            condition = true;
        end
    end
    % 目标信号积分旁瓣
    r = H_SL*x;
    Integrate_sidelobe(1,iter+1) = r' * r;
    % 干扰信号积分能量
    r_j = H_j*x;
    Integrate_jamming(1,iter+1) = r_j' * r_j;
    % 收敛条件
    current_total = Integrate_sidelobe(1,iter) + Integrate_jamming(1,iter);
    change_targetFunc = abs(current_total-prev_total);
    if change_targetFunc <= abs(current_total)*1e-4 
        Converg_Iter = iter;
    end
    prev_total = current_total;
    waitbar(iter/steps,a);
end
close(a);
if Converg_Iter ~= 0
    disp(['第',num2str(Converg_Iter),'次迭代收敛;']);
end
%% 设计优化波形后叠加干扰与噪声（时延）
signal_power = sum(abs(x).^2)/L;% 计算信号功率
interference_power = signal_power / (10^(SIR/20));% 计算干扰功率
noise_power = signal_power / (10^(SNR/20));% 计算噪声功率
% 噪声信号
noise = sqrt(noise_power/2)*(randn(L,1) + 1j * randn(L,1));
% 目标回波(多目标)
target1_echo = zeros(L,1);% 3200*1
target2_echo = zeros(L,1);% 3200*1
target1_delay_index = round(target1_delay * fs);
target2_delay_index = round(target2_delay * fs);
if target1_delay_index < length(t) && target2_delay_index < length(t)
    target1_echo(target1_delay_index:end,1) = RCS * x(1:end-target1_delay_index + 1,1);%目标1回波
    target2_echo(target2_delay_index:end,1) = RCS * x(1:end-target2_delay_index + 1,1);%目标2回波
end
% 干扰信号
x_j = x .* pulse_train;% 3200*1
Jammer_signal = x_j * sqrt(interference_power/sum(abs(x_j).^2)*L);
jammer_delay_index = round(jammer_delay * fs);%干扰机路程时延
%ISDRJ
jammer_ISDRJ = zeros(L,1);% 3200*1
Jammer_ISDRJ = zeros(L,1);% 3200*1
if jammer_delay_index < length(t)
    Jammer_ISDRJ(num_samps_jammer : end,1) = Jammer_signal(1 : end-num_samps_jammer+1,1);
    jammer_ISDRJ(jammer_delay_index:end,1) = Jammer_ISDRJ(1:end-jammer_delay_index + 1,1);
end
%ISPRJ
jammer_ISPRJ = zeros(L,1);% 3200*1
Jammer_ISPRJ = zeros(L,1);% 3200*1
for n = 0:num_periods-1% 第n次采样
    start_idx = round((n * Ts_jammer) * fs)+1;
    end_idx = start_idx + num_samps_jammer-1;
    Jammer_n = Jammer_signal(start_idx : end_idx,1);
    Jammer_ISPRJ(end_idx+1 : round(((n+1) * Ts_jammer) * fs),1) = [Jammer_n;Jammer_n;Jammer_n];
end
if jammer_delay_index < length(t)
    jammer_ISPRJ(jammer_delay_index:end,1) = Jammer_ISPRJ(1:end-jammer_delay_index + 1,1);
end
%ISCRJ
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
sig_single_ISDRJ = target1_echo + jammer_ISDRJ + noise;%最终回波3200*1
sig_single_ISPRJ = target1_echo + jammer_ISPRJ + noise;
sig_single_ISCRJ = target1_echo + jammer_ISCRJ + noise;
sig_multi_ISDRJ = target1_echo + target2_echo + jammer_ISDRJ + noise;
sig_multi_ISPRJ = target1_echo + target2_echo + jammer_ISPRJ + noise;
sig_multi_ISCRJ = target1_echo + target2_echo + jammer_ISCRJ + noise;
%% 仿真结果【单目标 2km】
h = conj(flip(h));
non_matched_output_single_ISDRJ = conv(sig_single_ISDRJ,h,'same');
r_xh = fftshift(non_matched_output_single_ISDRJ);
a = 1:length(r_xh);
b = (c*a)/(2*fs);
figure;
plot(b,abs(r_xh)/max(abs(r_xh)));
title('接收信号非匹配滤波输出');
xlabel('距离(m)');
ylabel('幅度(dB)');

non_matched_output_single_ISPRJ = conv(sig_single_ISPRJ,h,'same');
r_xh = fftshift(non_matched_output_single_ISPRJ);
a = 1:length(r_xh);
b = (c*a)/(2*fs);
figure;
plot(b,abs(r_xh)/max(abs(r_xh)));
title('接收信号非匹配滤波输出');
xlabel('距离(m)');
ylabel('幅度(dB)');

non_matched_output_single_ISCRJ = conv(sig_single_ISCRJ,h,'same');
r_xh = fftshift(non_matched_output_single_ISCRJ);
a = 1:length(r_xh);
b = (c*a)/(2*fs);
figure;
plot(b,abs(r_xh)/max(abs(r_xh)));
title('接收信号非匹配滤波输出');
xlabel('距离(m)');
ylabel('幅度(dB)');
%% 仿真结果【多目标 1.5km 2km】
non_matched_output_multi_ISDRJ = conv(sig_multi_ISDRJ,h,'same');
r_xh = fftshift(non_matched_output_multi_ISDRJ);
a = 1:length(r_xh);
b = (c*a)/(2*fs);
figure;
plot(b,abs(r_xh)/max(abs(r_xh)));
title('接收信号非匹配滤波输出');
xlabel('距离(m)');
ylabel('幅度(dB)');

non_matched_output_multi_ISPRJ = conv(sig_multi_ISPRJ,h,'same');
r_xh = fftshift(non_matched_output_multi_ISPRJ);
a = 1:length(r_xh);
b = (c*a)/(2*fs);
figure;
plot(b,abs(r_xh)/max(abs(r_xh)));
title('接收信号非匹配滤波输出');
xlabel('距离(m)');
ylabel('幅度(dB)');

non_matched_output_multi_ISCRJ = conv(sig_multi_ISCRJ,h,'same');
r_xh = fftshift(non_matched_output_multi_ISCRJ);
a = 1:length(r_xh);
b = (c*a)/(2*fs);
figure;
plot(b,abs(r_xh)/max(abs(r_xh)));
title('接收信号非匹配滤波输出');
xlabel('距离(m)');
ylabel('幅度(dB)');
%% 仿真结果
% 【实验1】验证可行
% 绘制设计的发射波形和非匹配滤波器时域波形【图2】
% z = 1:iter_max;
% figure;
% subplot(211);
% plot(z,real(Integrate_sidelobe(2:iter_max+1)), 'b-', 'LineWidth', 2);%旁瓣积分，蓝色实线，线宽为2
% hold on
% plot(z,real(Integrate_jamming(2:iter_max+1)), 'r--', 'LineWidth', 2);%干扰积分，红色虚线，线宽为2
% hold off
% title('目标函数收敛曲线');
% xlabel('迭代次数');
% ylabel('目标函数');
% legend('目标信号积分旁瓣', '干扰信号积分能量');
% subplot(212);
% plot(z,real(Integrate_sidelobe(2:iter_max+1) + Integrate_jamming(2:iter_max+1)));
% title('目标函数收敛曲线');
% xlabel('迭代次数');
% ylabel('目标函数');
% 
% figure;
% subplot(211);
% plot(z,20*log10(Integrate_sidelobe(2:iter_max+1)), 'b-', 'LineWidth', 2);%旁瓣积分，蓝色实线，线宽为2
% hold on
% plot(z,20*log10(Integrate_jamming(2:iter_max+1)), 'r--', 'LineWidth', 2);%干扰积分，红色虚线，线宽为2
% hold off
% title('目标函数收敛曲线(dB)');
% xlabel('迭代次数');
% ylabel('目标函数');
% legend('目标信号积分旁瓣', '干扰信号积分能量');
% subplot(212);
% plot(z,20*log10(Integrate_sidelobe(2:iter_max+1) + Integrate_jamming(2:iter_max+1)));
% title('目标函数收敛曲线(dB)');
% xlabel('迭代次数');
% ylabel('目标函数');

% % 绘制设计的发射波形和非匹配滤波器时域波形【图3】
% figure;
% subplot(2, 1, 1);
% plot(t*1e6,real(x));
% title('设计信号实部时域图');
% xlabel('时间(us)');
% ylabel('幅度(V)');
% subplot(2, 1, 2);
% plot(t*1e6,real(h));
% title('非匹配滤波器时域图');
% xlabel('时间(us)');
% ylabel('幅度(V)');

% % 计算并绘制信号和干扰信号的非匹配滤波输出【图4】
% h = conj(flip(h));
% non_matched_filter_output = conv(x,h);
% r_xh = 20 * log10(abs(non_matched_filter_output));
% figure;
% plot(2e6*t1,r_xh);
% title('发射信号非匹配滤波输出');
% xlabel('时间(us)');
% ylabel('幅度 (dB)');
% non_matched_filter_output_j = conv(x_j,h);
% r_xjh = 20 * log10(abs(non_matched_filter_output_j));
% figure;
% plot(2e6*t1,r_xjh);
% title('干扰信号非匹配滤波输出');
% xlabel('时间(us)');
% ylabel('幅度 (dB)');
% % 非匹配滤波峰值处理增益损耗
% x = x / norm(x);
% h = h / norm(h);
% LPG = 20 * log10((h'* x) / (h'*h));
% fprintf('非匹配滤波峰值处理增益损耗: %f dB\n', LPG);


