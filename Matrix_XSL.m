function [X_SL] = Matrix_XSL(L, LFM_signal)
    % 构建 X_SL 矩阵
    X_SL = zeros(L, 2*L-1);
    x = flip(LFM_signal');
    for i = 1:L  
        X_SL(i, i:L+i-1) = x;
    end
    X_SL(:,L) = 0;
    X_SL = X_SL';
end

