function [ind, x, d] = simplex_fase1_fase2(A, b, c, m, n)
    eps = 1e-9;

    for i = 1:m
        if b(i) < 0
            A(i, :) = -A(i, :);
            b(i) = -b(i);
        end
    end

    % ================================================================
    % FASE 1
    % ================================================================
    disp('=== Iniciando Fase 1 ===');

    A_aux = [A, eye(m)];
    c_aux = [zeros(n, 1); ones(m, 1)];

    ind_B = (n + 1):(n + m);
    ind_N = 1:n;

    [ind_fase1, x_aux, ind_B, ind_N] = ...
        executa_simplex(A_aux, b, c_aux, ...
                        ind_B, ind_N, m, n + m);
    
    % idx_fase1 = 0 // 
    % x_aux = solução otima


    val_fase1 = c_aux' * x_aux; % resposta otimal tem que ser 0
    % val_fase1 = 0 => x_aux_{idx_B} = 0

    if val_fase1 > eps                                          % usar 0 ou eps?
        disp('Resultado Fase 1: Problema inviavel.');

        ind = 1;
        x = [];
        d = [];
        return;
    end

    disp('Resultado Fase 1: Solucao viavel encontrada.');


    % ================================================================
    % TRANSICAO PARA A FASE 2
    % ================================================================
    disp('=== Iniciando Fase 2 ===');
    Aa = A; % Versao antiga de A para fazer a parte 2 do trabalaho
    % ^ TODO: Remove this later ^               % TODO: Remove this later <-

    if size(ind_B, 1) ~= m
        [A, ind_B, b, m] = pivotamento_neutro(A_aux, b, ind_B, m, n);
    end

    ind_B_fase2 = ind_B(ind_B <= n); % return only the elements <= n
    ind_N_fase2 = ind_N(ind_N <= n);

    B = A(:, ind_B_fase2);

    if size(B, 1) ~= m || size(B, 2) ~= m
        error(['ERRO CRITICO NA TRANSICAO PARA A FASE 2: ' ...
               'a matriz selecionada como base nao possui ' ...
               'dimensao m x m.']);
    end

    if rank(B) < m
        error(['ERRO CRITICO NA TRANSICAO PARA A FASE 2: ' ...
               'a matriz selecionada como base e singular.']);
    end


    % ================================================================
    % FASE 2
    % ================================================================
    [ind, x_full, ind_B_fase2, ind_N_fase2, d] = ...
        executa_simplex(A, b, c, ...
                        ind_B_fase2, ind_N_fase2, ...
                        m, n);

    if ind == -1
        disp('Resultado Fase 2: Problema ilimitado.');
    endif

    x = x_full;
    ind_B_fase2
    % AB = inv(Aa(:, ind_B_fase2))
    % cb = AB*c(ind_B_fase2)

    B_otima = Aa(:, ind_B_fase2);
    B_inv = inv(B_otima);

    % Vetor dual (y = (B^-1)^T * c_B)
    c_B = c(ind_B_fase2);
    y = B_inv' * c_B;

    % Vetor de custos reduzidos (c_bar = c - A^T * y)
    c_bar = c - Aa' * y;

    disp('Resultados da parte 2');
    disp('Inversa da Base Ótima (B^-1):');
    disp(B_inv);
    disp('Vetor Dual (y):');
    disp(y);
    disp('Vetor de Custos Reduzidos (c_bar):');
    disp(c_bar);

end


% ===================================================================
% SIMPLEX
% ===================================================================

function [ind, x_sol, ind_B, ind_N, d] = ...
    executa_simplex(A, b, c, ind_B, ind_N, m, n)

    eps = 1e-9;
    d = zeros(n, 1);

    max_iter = 100;
    iter = 0;

    while iter < max_iter

        iter = iter + 1;

        B = A(:, ind_B);
        N = A(:, ind_N);

        x_B = B \ b;
        % '\' operator means solving the Linear System Bx = b

        x_sol = zeros(n, 1);
        x_sol(ind_B) = x_B;

        y = (B') \ c(ind_B);

        c_N_bar = c(ind_N) - N' * y;

        % Teste de otimalidade
        if all(c_N_bar >= -eps)
            ind = 0;
            return;
        end

        % Variavel entrante
        cand_entra = find(c_N_bar <= -eps);

        idx_q = cand_entra(1);
        q = ind_N(idx_q);

        % Direcao simplex
        u = B \ A(:, q);

        % Teste de ilimitacao
        if all(u <= 0)
            ind = -1;

            d(q) = 1;
            d(ind_B) = -u;

            return;
        end

        % Teste da razao
        razoes = Inf(m, 1);

        for i = 1:m
            if u(i) > eps                                   % usar 0 ou eps?
                razoes(i) = x_B(i) / u(i);
            end
        end

        min_raz = min(razoes);

        cand_sai = find(razoes == min_raz);

        var_saindo = ind_B(cand_sai);

        [~, idx_s_local] = min(var_saindo);

        p_idx = cand_sai(idx_s_local);

        s = ind_B(p_idx);

        % Atualizacao da base
        ind_B(p_idx) = q;
        ind_N(idx_q) = s;

    end

    error('Numero maximo de iteracoes atingido.');

end
