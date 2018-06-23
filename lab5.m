clc; clear; close all;

function alpfa = ALPFA(max, n) 
        alpfa = log(1./max)./log(n);
end

function ret = Spline(x, y)
  x = x(:);
  n = length (x);
  a = y';

  b = c = zeros (size (a));
  h = diff (x);
  idx = ones (columns (a), 1);

  if (n == 2)
    b = (a(2,:) - a(1,:)) / (x(2) - x(1));
    a = a(1,:);
    d = [];
    c = [];
  elseif (n == 3)
    n = 2;
    c = (a(1,:) - a(3,:)) / ((x(3) - x(1)) * (x(2) - x(3))) ...
        + (a(2,:) - a(1,:)) / ((x(2) - x(1)) * (x(2) - x(3)));
    b = (a(2,:) - a(1,:)) * (x(3) - x(1)) ...
        / ((x(2) - x(1)) * (x(3) - x(2))) ...
        + (a(1,:) - a(3,:)) * (x(2) - x(1)) ...
        / ((x(3) - x(1)) * (x(3) - x(2)));
    a = a(1,:);
    d = [];
    x = [min(x), max(x)];
  else
    g = zeros (n-2, columns (a));
    g(1,:) = 3 / (h(1) + h(2)) ...
        * (a(3,:) - a(2,:) - h(2) / h(1) * (a(2,:) - a(1,:)));
    g(n-2,:) = 3 / (h(n-1) + h(n-2)) ...
        * (h(n-2) / h(n-1) * (a(n,:) - a(n-1,:)) - (a(n-1,:) - a(n-2,:)));

    if (n > 4)
      g(2:n - 3,:) = 3 * diff (a(3:n-1,:)) ./ h(3:n-2,idx) ...
          - 3 * diff (a(2:n-2,:)) ./ h(2:n - 3,idx);

      dg = 2 * (h(1:n-2) .+ h(2:n-1));
      dg(1) = dg(1) - h(1);
      dg(n-2) = dg(n-2) - h(n-1);

      ldg = udg = h(2:n-2);
      udg(1) = udg(1) - h(1);
      ldg(n - 3) = ldg(n-3) - h(n-1);
      c(2:n-1,:) = spdiags ([[ldg(:); 0], dg, [0; udg(:)]],
                            [-1, 0, 1], n-2, n-2) \ g;
    elseif (n == 4)
      dg = [h(1) + 2 * h(2); 2 * h(2) + h(3)];
      ldg = h(2) - h(3);
      udg = h(2) - h(1);
      c(2:n-1,:) = spdiags ([[ldg(:);0], dg, [0; udg(:)]],
                            [-1, 0, 1], n-2, n-2) \ g;
    endif

    c(1,:) = c(2,:) + h(1) / h(2) * (c(2,:) - c(3,:));
    c(n,:) = c(n-1,:) + h(n-1) / h(n-2) * (c(n-1,:) - c(n-2,:));
    b = diff (a) ./ h(1:n-1, idx) ...
        - h(1:n-1, idx) / 3 .* (c(2:n,:) + 2 * c(1:n-1,:));
    d = diff (c) ./ (3 * h(1:n-1, idx));

    d = d.'(:);
    c = c(1:n-1,:).'(:);
    b = b.'(:);
    a = a(1:n-1,:).'(:);
  endif

  ret = mkpp (x, cat (2, d, c, b, a));
end;


f = @(x) x.^2 .* sin(x);
df = @(x) 2*x .* sin(x) + x.^2 .* cos(x);
a = -1; b = 20;
ind = 1;

start = 200;
step = 200;
finish = 3000;

for N = start:step:finish
  x = linspace(a, b, N);
  y = f(x);
  h = x(2) - x(1);

  yf1 = (y(2:N) - y(1:N-1)) ./ h;
  yf2 = (y(3:N) - y(1:N-2)) ./ (2 * h);
  yf3 = (-3 .* y(1:N-2) + 4 .* y(2:N-1) - y(3:N)) ./ (2 * h);
  
  for i = 1:3:3*fix(N/3)
    t = x(i:i+2);
    p = Spline(t, f(t));
    p_diff = ppder(p);
    yfS(i:i+2) = ppval(p_diff, t);
  end;

  x1 = x(1:N-1);
  x2 = x(2:N-1);
  x3 = x(1:N-2);
  xS = x(1:3*fix(N/3));
  
  mod1(N/step) = max(abs(df(x1) - yf1));
  mod2(N/step) = max(abs(df(x2) - yf2));
  mod3(N/step) = max(abs(df(x3) - yf3));
  modS(N/step) = max(abs(df(xS) - yfS));
end;
  

figure;
hold on;
plot(x, y);
plot(x, df(x));
plot(x1, yf1);
plot(x2, yf2, 'k');
plot(x3, yf3, 'k');
plot(xS, yfS);
legend('original', 'original df', 'diff 1 near', 'diff 2 around', 'diff 3', 'diff 3 spline');

figure;
hold on;
plot(1:length(mod1), mod1)
plot(1:length(mod2), mod2)
plot(1:length(mod3), mod3)
plot(1:length(modS), modS)
legend('diff 1 near', 'diff 2 around', 'diff 3', 'diff 3 spline');

alpha1 = ALPFA(mod1, 1:length(mod1));
alpha2 = ALPFA(mod2, 1:length(mod2));
alpha3 = ALPFA(mod3, 1:length(mod3));
alphaS = ALPFA(modS, 1:length(modS));

figure;
hold on;
plot(1:length(mod1), alpha1);
plot(1:length(mod2), alpha2);
plot(1:length(mod3), alpha3);
plot(1:length(modS), alphaS);
legend('diff 1 near', 'diff 2 around', 'diff 3', 'diff 3 spline');

