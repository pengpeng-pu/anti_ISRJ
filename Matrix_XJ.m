function [X_j] = Matrix_XJ(L, pulse_train, LFM_signal)
    % 构建 X_j 矩阵
    X_j = zeros(L, 2*L-1);
    x = flip(LFM_signal');
    Z = flip(pulse_train') .* x;
    for i = 1:L  
        X_j(i, i:L+i-1) = Z;
    end
    X_j = X_j';
end








% 
% function [h] = mismatched_filter_design(x, x_j, Q, b_max, b_min)
%     % 提前计算 Q 的伪逆
%     Q_inv = pinv(Q + 1e-10 * eye(size(Q)));
%     % 计算一些中间项
%     x_Q_inv_x = x' * Q_inv * x;
%     x_j_Q_inv_x_j = x_j' * Q_inv * x_j;
%     x_Q_inv_x_j = x' * Q_inv * x_j;
%     % 计算分母
%     denominator = x_Q_inv_x * x_j_Q_inv_x_j - abs(x_Q_inv_x_j)^2;
%     if abs(denominator) < 1e-20
%         denominator = denominator + 1e-20;
%     end
%     % 计算分子部分
%     numerator1 = (b_max * x_j_Q_inv_x_j - b_min * x_Q_inv_x_j);
%     numerator2 = (b_min * x_Q_inv_x - b_max * x_Q_inv_x_j);
%     % 计算 h1 和 h2
%     h1 = (numerator1 / denominator) * Q_inv * x;
%     h2 = (numerator2 / denominator) * Q_inv * x_j;
%     % 计算最终结果
%     h = h1 + h2;
% end
% 
% function [x] = waveform_design(h_p, h_j_p, H_SL, H_j, x_prev, L, lambda, mu, nu, rho, b_max, b_min)
%     % 提前计算部分矩阵乘法
%     H_SLtH_SL = H_SL' * H_SL;
%     H_jtH_j = H_j' * H_j;
%     h_p_h_pT = h_p * h_p';
%     h_j_p_h_j_pT = h_j_p * h_j_p'; 
%     % 计算分子 1
%     numerator1 = 2 * H_SLtH_SL + 2 * H_jtH_j + rho * eye(L) + rho * h_p_h_pT + rho * h_j_p_h_j_pT;
%     % 为避免数值不稳定，添加一个小的正则化项
%     numerator1 = numerator1 + 1e-10 * eye(size(numerator1));
%     % 计算分子 2
%     numerator2 = rho * (x_prev - lambda) + rho * (b_max - mu) * h_p + rho * (b_min - nu) * h_j_p;
%     % 使用矩阵左除求解线性方程组
%     x = pinv(numerator1) * numerator2;
% end
% 
