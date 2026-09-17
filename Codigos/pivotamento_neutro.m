function [A_r,B_r,b, m] = pivotamento_neutro(A, b, B, m, n)
    inv_A = inv(A(:, B));
    A_can = inv_A * A;
    b = inv_A * b;

    remov = [];
    eps = 1e-9;

    for i = 1:size(B,2)
        if B(i) <= n
            continue;
        endif;

        cand = find(abs(A_can(i, :) - 0) >= eps);
        cand = cand(cand <= n);
        
        if (size(cand, 2) == 0)
            remov = [remov, i];
            m--;
            continue;
        endif
        
        pivo = min(cand);

        B_n = B;
        B_n(i) = pivo;

        inv_A = inv(A_can(:, B_n));
        A_can = inv_A * A_can;
        b = inv_A * b;
        B = B_n;
    end

    B_r = B;
    A_r = A_can;

    b(remov, :) = [];
    B_r(:, remov) = [];
    A_r(remov, :) = [];
    A_r = A_r(:, 1:n);

    return;
end;