%FMNF05 -- Project 1
%Authors:
%  Stefan Eng <atn08sen@student.lu.se>

%% Task 3

% Evaluate f(x) between -2 and 2 to get a feeling about the values.
xs = 0:0.1:1;
ys = zeros(size(xs));

for i = 1:numel(xs)
  ys(i) = f(xs(i));
end

fplot('t3_check', xs, ys, '', 'f(x)', 10, 3);

% Test the bisection method.
[xa, iters] = bisection(@f, 0, 1, 0.5*10^(-3), 100);
fxa = f(xa);
fprintf('Task3 -- Bisection\n');
fprintf('---\n');
fprintf('Bisection method produces: f(%.5f) = %f in %d iterations.\n', xa, fxa, iters);

% Test the fixed point method for cos.
%[xa, iters] = fixed_point(@cos, 0, 0.5*10^-9, 100)
k = xa;

% Run fixed-point for g(x), with xa from prev as guess.
fprintf('\n');
fprintf('Task4 -- Fixed point.\n');
fprintf('---\n');
[xa, iters] = fixed_point(@g2, xa, 0.5*10^-6, 100);
fprintf('fixed-point produces xa: %.6f in %d iterations.\n', xa, iters);
fprintf('with value f(%.6f): %.16f.\n', xa, f(xa));
kc = xa;

%exp(-k)*(20+k);
%fprintf('Value of g2''(x): %.6f\n', 6/(k^2 + 34*k + 128));

fprintf('\n');
fprintf('Task5 -- Newton-Raphson.\n');
fprintf('---\n');
eps_newton = 0.5*10^-17;
[xa, iters] = newtonraphs(@ft, @ftp, kc, -1.25, eps_newton, 500);
fprintf('Newton-Raphson produces xa: %.16f in %d iterations.\n', xa, iters);
fprintf('with value f(%.6f): %.16f.\n', xa, ft(xa, kc));

xs = -10:0.1:10;
ys = zeros(size(xs));
for i = 1:numel(xs)
  ys(i) = ft(xs(i), kc);
end
fplot('t5_check', xs, ys, '', 'f(t)', 10, 3);

xs = -2:0.1:0;
ys = zeros(size(xs));
for i = 1:numel(xs)
  ys(i) = ft(xs(i), kc);
end
fplot('t5_check_2', xs, ys, '', 'f(t)', 10, 3);

fprintf('\n');
fprintf('Task6 -- Compute tc with bisection and fixed-point.\n');
fprintf('---\n');

ftssolve =@(x) ft(x, kc);
[xa, iters] = bisection(ftssolve, -1.5, 0, eps_newton, 2600);
fprintf('Bisection method computed xa = %.16f in %d iterations.\n', xa, iters);
fprintf('Bisection error: %.16f.\n', ft(xa, kc));

fixed_g3 =@(x) g3(x, kc);
[xa, iters] = fixed_point(fixed_g3, -1.2, eps_newton/10, 5000);
fprintf('Fixed-point method computed xa = %.16f in %d iterations.\n', xa, iters);
fprintf('Fixed-point error: %.16f.\n', ft(xa, kc));

%g3d(-1.2, kc)

%xs = -2:0.1:0;
%ys = zeros(size(xs));
%for i = 1:numel(xs)
%  ys(i) = g3(xs(i), kc);
%end
%fplot('t6_check_1', xs, ys, '', 'f(t)', 10, 3);


fprintf('\n');
fprintf('Task7 -- Comapare to fzero results.\n');
fprintf('---\n');

fkc = fzero(@f, 1);
fprintf('fzero solution for \tkc = %.16f.\n', fkc);
fprintf('our solution for \tkc = %.16f.\n', kc);
tc = fzero(ftssolve, 1);
fprintf('fzero solution for tc with our kc, \t\ttc = %.16f.\n', tc);
ftssolve =@(x) ft(x, fkc);
tc = fzero(ftssolve, -1.1);
fprintf('fzero solution for tc with fsolve kc, \t\ttc = %.16f.\n', tc);
fprintf('Fixed-point error: %.16f.\n', ft(t, fkc));
fprintf('T(t) Result for tc = %.16f.\n', ft(tc, kc));

function y = f(x)
  y = 7 + 0.5 * x - (10 + 0.5 * x ) * exp(-x);
end

function y = g(x)
  y = -2*(10 + 0.5*x)*exp(-x)+14;
end

function y = g2(x)
  y = log((-10 - 0.5*x)/(-7-0.5*x));
end

function y = ft(t, k)
  y = -15 + 0.5*t - 0.5*k + (10 + 0.5*k)*exp(-k*t);
end

function y = ftp(t, k)
  y = 0.5 - k*(10 + 0.5*k)*exp(-k*t);
end

function y = g3(t, k)
  y = log((15-0.5*t+0.5*k)/(10+0.5*k)) * (1/-k);
end

function y = g3d(t, k)
  y = -1/(k*(k-t+30));
end

% error = | xa - r | < (b-a)/2^(n+1)
% Correct to p decimal places -> error < 0.5 * 10^(-p)

function [xc, iters] = bisection(f, a, b, eps, Nmax)
  % Calculate function values based on inverval.
  fa = f(a); fb = f(b);
  iters = 1;
  while iters <= Nmax
    % Break if eps is satisfied.
    if (b-a)/2 < eps; break; end
    % Caculate new mid and function value.
    mid = (a+b)/2; fmid = f(mid);
    % Break if fmid is a root.
    if fmid == 0; break; end
    % Not root, construct new interval.
    if fa * fmid < 0 % [a, mid] includes root.
      b = mid;
    else % [mid, b] includes root.
      a = mid; fa = fmid;
    end
    iters = iters + 1;
  end
  % Calculate and return approximate root.
  xc = (a+b)/2;
end

% Construct fixed point algorithm.

function [xc, iters] = fixed_point(g, guess, eps, Nmax)
  iters = 1;
  xc = guess;
  while iters <= Nmax
    % Generate next guess based on current guess.
    nxc = g(xc);
    % Calculate difference between previous and current guess.
%    diff = abs(((nxc - xc)*(nxc + xc))/(nxc + xc));
    diff = abs(xc - nxc);
    % Use current nxc value as next xc value.
    fprintf('%.19f %.19f %.19f %.25f\n', xc, nxc, diff, eps);
    xc = nxc;
    % Break if new guess similar enough to previous guess.
    if diff < eps; break; end
    iters = iters + 1;
  end
end

function [xc, iters] = newtonraphs(f, fp, k, guess, eps, Nmax)
  iters = 1;
  xc = guess;
  while iters <= Nmax
    nxc = xc - (f(xc, k)/fp(xc, k));
    err = f(nxc, k);
    iters = iters + 1;
    xc = nxc;
    if abs(err) < eps; break; end
  end
end


function fplot(name, xs, ys, x_label, y_label, width, height)
  fig = figure('visible','off');
  plot(xs, ys);
  set(gcf,'Units','centimeters');
  screenposition = get(gcf,'Position');
  set(gcf,...
      'PaperPosition',[0 0 width height],...
      'PaperSize', [width height]);
  set(gca, 'FontSize', 8);
  set(gca, 'FontName', 'Computer Modern');
  xlabel(x_label);
  ylabel(y_label);
  saveas(fig, ['figs/', name], 'pdf');
end
