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

function [P, yy_P] = REMEZ(f, x)
  N = l = length(x);
  E = 1e-6;
  cur_ksi = x;
  xx = linspace(x(1), x(end), 100);
  yy = f(xx);
  
  while (1)
    y = f(x);
    sl_r = zeros(N);
    sl_l = y';

    for i=1:N
      sl_r(i,1) = 1;
      
      for(j=1:N-2)
        sl_r(i, j + 1) = x(i)^j;
      end;
      
      sl_r(i, N) = (-1)^(i - 1);
    end;

    P = fliplr((sl_r\sl_l)');
    d = abs(P(N));
    P = P(2:end);
    yy_P = polyval(P, xx);
    [D, Di] = max(abs(yy - yy_P));
    cycle = 0;
    
    for c_i=1:length(cur_ksi)
      if (xx(Di) ==  cur_ksi(c_i))
        cycle = 1;
        break;
      end;
    end;
    
    if(cycle == 1)
      break;
    end;
    
    cur_ksi(l) = xx(Di);

    for i=1:N
      x_m(i) = abs(xx(Di) - x(i));
    end;
    
    [x_mv, x_mi]  = min(x_m);
    
    if (x_mi != 1 && x_mi != N)
        x(x_mi) = xx(Di);
    end;
    
    if (x_mi == N)
       new_x(1) = x(1);
       new_x(N) = x(N);
       
       for j = 2:N-2
         new_x(j) = x(j + 1);
       end;
       
       new_x(N - 1) = xx(Di);
       x = new_x;
    end;
      if (x_mi == 1)
       new_x(1) = x(1);
       new_x(N) = x(N);
       
       for j = 2:N-2
         new_x(j + 1) = x(j);
       end;
       
       new_x(2) = xx(Di);
       x = new_x;
    end;
    
    l = l + 1;
    
    if (D - d < E )
      break;
    end;
  end;
end;

f = @(x)  10 * e.^-x;
F = @(x)  - 10 * e.^-x;

a = 0; b = 5;

Rdots = 5;
Sdots = 3;

left = 201; right = 1001;
step = 200;

true_integr = F(b) - F(a);

count = 1;

for N = left:step:right
  x = linspace(a, b, N);
  y = f(x);
  delta = x(2) - x(1);

  integr1 = delta .* sum(y(2:N));  
  integr2 = (delta/2) .* sum(y(2:N) + y(1:N-1));
  
  true_integrSimp = F(b) - F(a);
  if mod(N, 2) == 1
    integr3 = (delta / 3) * (y(1) + 4 * sum(f(x(2:2:end-1))) + ...
                                 2 * sum(f(x(3:2:end-1))) + y(end));
  else
    integr3 = (delta / 3) * (y(1) + 4 * sum(y(2:2:end-1)) + ...
                                 2 * sum(y(3:2:end-2)) + y(end-1));
                                 
      true_integrSimp = F(x(end-1)) - F(a);
  end;
                        
  
  integrR = 0;
  for i = 1:Rdots-1:N
    if (i+Rdots-1 > N)
      break;
    end;
    t = x(i:i+Rdots-1);
    p = REMEZ(f, t);
    pp = polyint(p);
    integrR += (polyval(pp, t(end)) - polyval(pp, t(1)));
  end;
  true_integrR = F(t(end)) - F(a);
  
  integrS = 0;
  for i = 1:Sdots-1:N
    if (i+Sdots-1 > N)
      break;
    end;
    t = x(i:i+Sdots-1);
    p = Spline(t, f(t));  
    pp = ppint(p);
    integrS += (ppval(pp, t(end)) - ppval(pp, t(1)));
  end
  true_integrS = F(t(end)) - F(a);
  
  mod1(count) = abs(true_integr - integr1);
  mod2(count) = abs(true_integr - integr2);
  mod3(count) = abs(true_integrSimp - integr3);
  modR(count) = abs(true_integrR - integrR);
  modS(count) = abs(true_integrS - integrS);
  
  count++;
end;


disp("F(b) - F(a) =  "); disp(true_integr)
disp("1 dot: "); disp(integr1)
disp("2 dots: "); disp(integr2)
disp("3 dots(Simpson): "); disp(integr3)
disp("Remez: ");disp(integrR)
disp("Spline: ");disp(integrS)

dist = left:step:right;

figure;
hold on;
semilogy(dist, mod1, 'r');
semilogy(dist, mod2, 'k');
semilogy(dist, mod3, 'm');
semilogy(dist, modR, 'g');
semilogy(dist, modS, 'b');
legend('1 dot', '2 dots', '3 dots (Simpson)', 'Remez', 'Spline');


alpha1 = ALPFA(mod1, 1:length(mod1));
alpha2 = ALPFA(mod2, 1:length(mod2));
alpha3 = ALPFA(mod3, 1:length(mod3));
alphaR = ALPFA(modR, 1:length(modR));
alphaS = ALPFA(modS, 1:length(modS));

figure;
hold on;
plot(dist, alpha1, 'r');
plot(dist, alpha2, 'k');
plot(dist, alpha3, 'm');
plot(dist, alphaR, 'g');
plot(dist, alphaS, 'b');
legend('1 dot', '2 dots', '3 dots (Simpson)', 'Remez', 'Spline');

