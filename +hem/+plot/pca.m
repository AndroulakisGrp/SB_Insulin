function pca(t_plot, y_plot, vars, pc, mu, sigma)
% Plot multiple curves in PC space.

M = size(y_plot.P, 1); % Number of profiles

hold on;
% Plot each individual
for i=1:M
% for i=[2, 22, 24]
    z = zeros(length(t_plot), length(vars));
    % z(k,j) = values for one simulation, the jth variable at the kth time point
    for j=1:length(vars)
        z(:,j) = eval(sprintf('y_plot.%s(i,:)', vars{j}))';
    end
%     zz = zscore(z);
    zz = bsxfun(@minus, z, mu);
    zz = bsxfun(@rdivide, zz, sigma);
    
    % Phase plots with lines colored by individual
%     plot(zz*pc(:,1), zz*pc(:,2), 'Color', colors(i,:));

    % Phase plots with lines colored by eigenvalue
    x = zz*pc(:,1);
    y = zz*pc(:,2);
    % Find points above the level
    level1 = (y_plot.e_d(i,:)>=-1000);
    level2 = (y_plot.e_d(i,:)>=0.02);
    level3 = (y_plot.e_d(i,:)>=0.3);
    level4 = (y_plot.e_d(i,:)>=1.5);
    % Create 4 copies of y
    line1 = y;
    line2 = y;
    line3 = y;
    line4 = y;
    % Set the values you don't want to get drawn to nan
    line1(~level1) = NaN;
    line2(~level2) = NaN;
    line3(~level3) = NaN;
    line4(~level4) = NaN;
    plot(x, line1, 'g', x, line2, 'y', x, line3, 'r', x, line4, 'k');
end
xlabel('PC1');
ylabel('PC2');