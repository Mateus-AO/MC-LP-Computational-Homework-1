clc; clear;

d1 = 3;
d2 = 3;
d3 = 5;
d4 = 0;

A = [2, 1, 1;
    1, 2, d4+1;
    1, 1, 2];

A = [A, eye(3)];

b = [ 10+d1;
      12+d2;
      8+d3
    ];

c = [-(d1+1); -(d2+2); -(d3+1); 0; 0; 0];

[m, n] = size(A);

save('instancia_matricula.mat', 'A', 'b', 'c', 'm', 'n');