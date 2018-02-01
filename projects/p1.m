%FMNF05 -- Project 1
%Authors:
%  Stefan Eng <atn08sen@student.lu.se>

%% Task 3

xs = -2:0.1:2;
ys = zeros(size(xs));

for i = 1:numel(xs)
  ys(i) = f(xs(i));
end

f(-2)
f(2)

fplot('t3_check', xs, ys, '', 'f(x)', 4, 3);

ys;
function y = f(x)
y = 7 + 0.5 * x - (10 + 0.5 * x ) * exp(-x);
end

function fplot(name, xs, ys, x_label, y_label, width, height)
  fig = figure('visible','off');
  plot(xs, ys);
  set(gcf,'Units','centimeters');
  screenposition = get(gcf,'Position')
  set(gcf,...
      'PaperPosition',[0 0 width height],...
      'PaperSize', [width height]);
  set(gca, 'FontSize', 8);
  set(gca, 'FontName', 'Computer Modern');
  xlabel(x_label);
  ylabel(y_label);
  saveas(fig, ['figs/', name], 'pdf');
end
