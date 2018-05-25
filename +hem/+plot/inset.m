function inset(t, y, pos, hoffset, color)
% Plots (t, y) on new axes at pos, along with a vertical bar at hoffset.
% NOTE: This function should probably be generalized, particularly the
% rectangle part

h3 = axes('Position', pos);
hem.plot.var(t, y, [], [], 48, [], 'k');
hold on;
% plot(48+[hoffset, hoffset], [0, 3], 'b', 'LineWidth', 6);
rectangle('Position', [hoffset, 0, 1, 3], 'FaceColor', color, 'LineStyle', 'none')
set(gca, 'XTick', []);
set(gca, 'YTick', []);
box on;
axis([0, 24, 0, 3]);