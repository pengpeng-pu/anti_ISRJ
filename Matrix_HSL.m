function [H_SL] = Matrix_HSL(L, h)
    % 构建 H_SL 矩阵
    H_SL = zeros(L, 2*L-1);
    h_trans = h';
    for i = 1:L  
        H_SL(i, L+1-i:L+L-i) = h_trans;
    end
    H_SL(:,L) = 0;
    H_SL = H_SL';
end
