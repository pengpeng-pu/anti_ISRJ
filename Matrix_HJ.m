function [H_j] = Matrix_HJ(L, pulse_train, h)
    % 构建 H_j 矩阵
    H_j = zeros(L, 2*L-1);
    h_trans = h';
    for i = 1:L  
        H_j(i, L+1-i:L+L-i) = pulse_train(i) * h_trans;
    end
    H_j = H_j';
end
