function var(t, y, data_t, data_y, t_pre, line_style, color, varargin)
%HEM.PLOT.VAR Plot one variable and data.
%   This function is usually called from HEM.PLOT.VARS rather than
%   directly.
%
%   HEM.PLOT.VAR(T,Y) plots vector Y versus vector T.
%
%   HEM.PLOT.VAR(T,Y,DATA_T,DATA_Y) also plots discrete data points defined
%   in the vectors DATA_T and DATA_Y.
%
%   HEM.PLOT.VAR(T,Y,DATA_T,DATA_Y,T_PRE) delays the plotting of the line
%   and the data points by T_PRE hours.
%
%   HEM.PLOT.VAR(T,Y,DATA_T,DATA_Y,T_PRE,LINE_STYLE) plots the line (T,Y)
%   withthe line type defined in LINE_STYLE; for instance, '--' for a
%   dashed line (see PLOT).
%
%   HEM.PLOT.VAR(T,Y,DATA_T,DATA_Y,T_PRE,LINE_STYLE,COLOR) plots both the
%   line and the data points with the color set in COLOR ('r', 'g', 'b',
%   etc.; see PLOT).
%
%   Any further arguments will be passed on to PLOT for the plot of Y
%   versus T.
%
%   See also HEM.PLOT.VARS, PLOT.

% Default values
if nargin<=2
    data_t = [];
end
if nargin<=3
    data_y = [];
end
if nargin<=4
    t_pre = 48;
end
if nargin<=5
    line_style = '';
end
if nargin<=6
    color = 'k';
end

% Error checking
if length(data_t) ~= length(data_y)
    error('data_t and data_y must be vectors of equal length.');
end
if length(t) ~= length(y)
    error('t and y must be vectors of equal length.');
end

% Plot stuff
hold on;
plot(t(t>t_pre)-t_pre, y(t>t_pre), sprintf('%s%s', color, line_style), 'LineWidth', 1, varargin{:})
if ~isempty(data_t) && ~isempty(data_y)
    plot(data_t, data_y, sprintf('%s.', color), 'MarkerSize', 20);
end
hold off;

% y-axis scaling
% ys = [data_y, y(t>t_pre)'];
% y0 = min(ys);
% y1 = max(ys);
% ylim([y0-0.1*(y1-y0), y1+0.1*(y1-y0)]);
