function [x] = waveform_design(x_numerator1, h_p, h_j_p, x_prev, lambda, mu, nu, rho, b_max, b_min)
%     % 提前计算部分矩阵乘法
%     H_SLtH_SL = H_SL' * H_SL;
%     H_jtH_j = H_j' * H_j;
%     h_p_h_pT = h_p * h_p';
%     h_j_p_h_j_pT = h_j_p * h_j_p'; 
%     % 计算分子 1
%     numerator1 = 2 * H_SLtH_SL + 2 * H_jtH_j + rho * eye(L) + rho * h_p_h_pT + rho * h_j_p_h_j_pT;
%     % 为避免数值不稳定，添加一个小的正则化项
%     numerator1 = numerator1 + 1e-10 * eye(size(numerator1));
    % 计算分子 2
    numerator2 = rho * (x_prev - lambda) + rho * (b_max - mu) * h_p + rho * (b_min - nu) * h_j_p;
    % 使用矩阵左除求解线性方程组
    x = x_numerator1 * numerator2;
end