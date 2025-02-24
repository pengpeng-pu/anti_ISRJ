function [h] = mismatched_filter_design(x, x_j, Q_inv, b_max, b_min)
    % 计算一些中间项
    x_Q_inv_x = x' * Q_inv * x; %1
    x_j_Q_inv_x_j = x_j' * Q_inv * x_j; %
    x_Q_inv_x_j = x' * Q_inv * x_j;
    % 计算分母
    denominator = x_Q_inv_x * x_j_Q_inv_x_j - abs(x_Q_inv_x_j)^2;
    if abs(denominator) < 1e-20
        denominator = denominator + 1e-20;
    end
    % 计算分子部分
    numerator1 = (b_max * x_j_Q_inv_x_j - b_min * x_Q_inv_x_j);
    numerator2 = (b_min * x_Q_inv_x - b_max * x_Q_inv_x_j);
    % 计算 h1 和 h2
    h1 = (numerator1 / denominator) * Q_inv * x;
    h2 = (numerator2 / denominator) * Q_inv * x_j;
    % 计算最终结果
    h = h1 + h2;
end