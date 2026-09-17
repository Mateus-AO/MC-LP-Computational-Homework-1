% GERAR_INSTANCIA_CRITICA.M
% Script para gerar a instância de teste onde o Simplex ingênuo falha.

clc; clear;

% Problema com restrição redundante para forçar x_artifical = 0 na base final da Fase 1
A = [ 1,  2,  0;
      2,  4,  0;
      1,  1,  1 ];

b = [ 4;
      8;
      3 ];

c = [-1; -2; 0];

[m, n] = size(A);

save('instancia_critica.mat', 'A', 'b', 'c', 'm', 'n');

disp('Arquivo "instancia_critica.mat" gerado com sucesso!');
disp('Para testar a falha no Octave execute:');
disp('  load("instancia_critica.mat");');
disp('  [ind, x, d] = simplex_fase1_fase2(A, b, c, m, n);');
