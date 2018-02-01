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
[xa, iters] = bisection(@f, 0, 1, 0.5*10^(-2), 100);
fxa = f(xa);
fprintf('Bisection method produces: f(%.5f) = %f in %d iterations.\n', xa, fxa, iters);

ys;
function y = f(x)
y = 7 + 0.5 * x - (10 + 0.5 * x ) * exp(-x);
end

% error = | xa - r | < (b-a)/2^(n+1)
% Correct to p decimal places -> error < 0.5 * 10^(-p)

function [xc, iters] = bisection(f, a, b, eps, Nmax)
  % Calculate function values based on inverval.
  fa = f(a); fb = f(b);
  iters = 1;
  while iters <= Nmax
    % Break if eps is satisfied.
    if (b-a) < eps; break; end
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
