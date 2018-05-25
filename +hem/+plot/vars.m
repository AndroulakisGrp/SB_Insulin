function vars(t, y, varargin)
%HEM.PLOT.VARS Plot multiple model variables.
%   HEM.PLOT.VARS(T,Y) plots some fields of the structure Y in subplots as
%   functions of the vector T. T and Y are typically generated as the
%   output of HEM.MODEL.RUN. Given no further arguments, F, EPI, P, and A
%   will be plotted in a 2 by 2 subplot.
%
%   The arguments T,Y can be followed by parameter/value paris to specify
%   additional options (these can also be passed in a struct):
%
%      'vars_to_plot' is a cell array containing variable names to plot
%      'all_vars' is a boolean. When true, it overrides var_to_plot and
%         plots all variables in y
%      'titles' is a cell array containing titles of subplots
%      't_pre' is the number of hours elapsed before plotting begins
%         (should be same as passed to HEM.MODEL.RUN)
%      'dim' is a vector containing the dimensions of the subplot [nx,ny]
%      'line_style' is a string containing the line type for the curves;
%         for instance, '--' for a dashed line (see PLOT)
%      'color' is a string containing the color of the curves (see PLOT)
%      'plot_data' is a string set to the name of a file in the +data
%         folder to plot data along with the model output
%      'plot_data_offset' is the hour at which to start plotting data
%
%   See also HEM.PLOT.VAR, HEM.MODEL.RUN, PLOT.

%% Function parameters
p = inputParser;
p.StructExpand = true;

% Which variables to plot
p.addParamValue('vars_to_plot', {'F', 'EPI', 'P', 'A'}, @(x)validateattributes(x, {'cell'}, {}));
% Or... override the previous line and plot all variables?
p.addParamValue('all_vars', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% Titles of the variables to plot
p.addParamValue('titles', {}, @(x)validateattributes(x, {'cell'}, {}));
% How many hours elapse before plotting begins?
p.addParamValue('t_pre', 48, @(x)validateattributes(x, {'numeric'}, {'scalar'}));
% Dimension of subplot
p.addParamValue('dim', [0, 0], @(x)validateattributes(x, {'numeric'}, {'vector', 'size', [1, 2]}));
% Line style for the main plot
p.addParamValue('line_style', '-', @(x)validateattributes(x, {'char'}, {}));
% Color for the main plot
p.addParamValue('color', 'k', @(x)validateattributes(x, {'char'}, {}));
% To plot data, set to the name of the file in the +data folder
p.addParamValue('plot_data', '', @(x)validateattributes(x, {'char'}, {}));
% What hour to start plotting data at
p.addParamValue('plot_data_offset', 0, @(x)validateattributes(x, {'numeric'}, {'scalar'}));

p.parse(varargin{:});
s = p.Results;

if s.all_vars
    s.vars_to_plot = fieldnames(y)';
end

if isempty(s.titles) || s.all_vars
    s.titles = s.vars_to_plot;
end

% BMES stochastic poster
% if y.P(end) > 5
%     s.color = 'r';
% end

n = length(s.vars_to_plot);

% Make appropriately sized subplot
if sum(s.dim) == 0
    nx = ceil(sqrt(n));
    ny = ceil(sqrt(n));
    if nx*(ny-1) >= n
        ny = ny - 1;
    end 
else
    nx = s.dim(1);
    ny = s.dim(2);
end

for i=1:n
    subplot(nx, ny, i);
    if ~isempty(s.plot_data)
        [data_t, data_y] = hem.data.(s.plot_data)(s.vars_to_plot{i});
        data_t = data_t + s.plot_data_offset;
%         % Normalize for starting time - NOTE: This is kind of BS
%         if ~isempty(data_t)
%             data_y = data_y*y.(vars_to_plot{i})(find(t>data_t(1), 1));
%         end
    else
        data_t = [];
        data_y = [];
    end
    if strcmp('Fboth', s.vars_to_plot{i})
        y_plot = y.F + y.F_ex;
    else
        y_plot = y.(s.vars_to_plot{i});
    end

    hem.plot.var(t, y_plot, data_t, data_y, s.t_pre, s.line_style, s.color);
    ylabel(s.titles{i});

%     set(gca, 'XTick', [0:12:48]);
%     set(gca, 'XTickLabel', {'12am', '12pm'})
%     xlim([0, 48]);
    grid on;
%     suplabel('Time (hr)', 'x', [.1 .1 .85 .85]);
end

% if s.extra_plots
%     % Plot all variables along with their max and min values
%     figure;
%     vars = {'LPS', 'EPI', 'F', 'HRV', 'M', 'mRNA_R', 'R', 'LPSR', 'IKK', 'mRNA_IkBa', 'IkBa', 'NFkBn', 'P', 'A', 'E', 'R_EPI', 'EPIR', 'cAMP', 'mRNA_R_F', 'R_F', 'FR', 'FRn', 'TC', 'A1', 'A2', 'T_sym', 'T_par', 'F_ex', 'mRNA_R_F_en', 'R_F_en', 'FR_en', 'FRn_en'};
%     n = length(vars);
%     for i=1:n
%         var = vars{i};
%         subplot(ceil(sqrt(n)), ceil(sqrt(n)), i);
%         plot(t, mean(y.(var), 1));
%         hold on;
%         title(var);
%     end
% end
