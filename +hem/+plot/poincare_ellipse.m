function [SD1, SD2] = poincare_ellipse(x0, c)
%HEM.PLOT.POINCARE_ELLIPSE Plot an ellipse representing a poincare map.
%   HEM.PLOT.POINCARE_ELLIPSE(X) plots an ellipse where the axis radiuses
%   are the standard deviations of X along the diagonal lines y=x and y=-x.
%   Based on (for example) "Do Existing Measures of Poincaré Plot Geometry
%   Reflect Nonlinear Features of Heart Rate Variability?" by Brennan, et
%   al.
%
%   HEM.PLOT.POINCARE_ELLIPSE(X, C) makes the ellipse have the color
%   defined by the RGB vector in C.
%
%   [SD1,SD2] = HEM.PLOT.POINCARE_ELLIPSE(X) returns the standard
%   deviations along the diagonal axes that represent the ellipse radiuses.
%   Similar to MATLAB's HIST function, if you ask for these SD1 and SD2
%   return values, the ellipse is not plotted.
%
%   See also HEM.MODEL.IPFM, HEM.PLOT.POINCARE.

if size(x0,1) ~= 1
    x0 = x0';
end
if size(x0,1) ~= 1
    error('X must be a vector.');
end

if nargin < 2
    c = [0.8, 0.8, 0.8];
end

theta = pi/4; % Rotate 45 degrees

x = [cos(theta), -sin(theta); sin(theta), cos(theta)]*[x0(1:end-1); x0(2:end)];

SD1 = std(x(1,:));
SD2 = std(x(2,:));

if nargout == 0
    % figure;
    % plot_poincare(x);
    
    % http://www.mathworks.com/matlabcentral/fileexchange/289-ellipse-m
    cx = mean(x0(1:end-1));
    cy = mean(x0(2:end));
    h = ellipse(SD2, SD1, theta, cx, cy);
    
    % Draw radii
    h2 = plot([cx-SD1/sqrt(2), cx+SD1/sqrt(2)], [cy+SD1/sqrt(2), cx-SD1/sqrt(2)], 'LineWidth', 2);
    h3 = plot([cx-SD2/sqrt(2), cx+SD2/sqrt(2)], [cy-SD2/sqrt(2), cx+SD2/sqrt(2)], 'LineWidth', 2);

    set(h, 'LineWidth', 3);
    set(h, 'Color', c);
    set(h2, 'Color', c);
    set(h3, 'Color', c);
    
    clear SD1 SD2; % Suppress outpt
end

% title(sprintf('%g %g', SD1, SD2));

% figure;
% subplot(1,2,1);
% hist(x(1,:), 100);
% subplot(1,2,2);
% hist(x(2,:), 100);