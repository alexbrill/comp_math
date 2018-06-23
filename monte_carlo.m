clc; clear; close all;

rand("seed", 1);

f = @(x)  3 * x .^ 2 + 2 * x + 1;
F = @(x)  x .^ 3 + x .^ 2 + x;

a = 0; b = 0.5;

FF = F(b) - F(a);

N_s = 1000; s = 1000; N_f = 200000;
for N = N_s:s:N_f
  X = a + (b - a) * rand(1, N);  
  F1(N/s) = (b - a) * sum(f(X)) / N;
  mod(N/s) = abs(FF - F1(N/s));
  alpha(N/s) = log(1/mod(N/s)) / log(N);
end

disp("F(b) - F(a) =  "); disp(FF);
[min_mod, ind] = min(mod);
disp("Best result: "); disp(F1(ind));
disp("with mod: ");disp(min_mod)

figure;
plot(N_s:s:N_f, F1);

figure
semilogy(N_s:s:N_f, mod, 'LineWidth', 2)
title('mod');
grid on

figure
plot(N_s:s:N_f, alpha, 'LineWidth', 2)
title('alpha');
grid on
