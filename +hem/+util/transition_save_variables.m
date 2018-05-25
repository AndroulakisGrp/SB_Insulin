function [t_plot, y_plot] = transition_save_variables(LPS_0s, t_LPSs, varargin)
% This is basically a convenience function to do a bunch of crap that is
% needed in HEM.SCRIPTS.TRANSITION, namely running a bunch of simulations
% while varying either the initial LPS dose or the time of LPS exposure,
% set by passing a scalar as one of the first two arguments and a vector as
% the other. There are also a bunch of options that can be set.

p = inputParser;
addRequired(p, 'LPS_0s');
addRequired(p, 't_LPSs');
addParamValue(p, 'eig', true, @islogical);
addParamValue(p, 'init_parfor_progress', true, @islogical);
addParamValue(p, 't_ex', [0, 0], @(x)validateattributes(x, {'numeric'}, {'vector', 'size', [1, 2]}));
addParamValue(p, 't_plot', linspace(-24, 48, 1000), @(x)validateattributes(x, {'numeric'}, {'vector'}));
addParamValue(p, 'center_plots', true, @islogical);
parse(p, LPS_0s, t_LPSs, varargin{:});
eig = p.Results.eig;
init_parfor_progress = p.Results.init_parfor_progress;
t_ex = p.Results.t_ex;
t_plot = p.Results.t_plot;
center_plots = p.Results.center_plots;

if length(LPS_0s) == 1 && length(t_LPSs) > 1
    vary = 't_LPS';
    M = length(t_LPSs);
elseif length(LPS_0s) > 1 && length(t_LPSs) == 1
    vary = 'LPS_0';
    M = length(LPS_0s);
elseif length(LPS_0s) == 1 && length(t_LPSs) == 1
    warning('Both LPS_0s and t_LPSs (the first two arguments) are scalars, so only one simulation will be run.');
    vary = 't_LPS';
    M = 1;
else
    error('Bad inputs.');
end

vars = {'LPS', 'LPSR', 'IKK', 'mRNA_R', 'R', 'NFkBn', 'mRNA_IkBa', 'IkBa', 'P', 'A', 'E', 'EPI', 'cAMP', 'M', 'F', 'mRNA_R_F', 'R_F', 'FR', 'FRn', 'R_EPI', 'EPIR'};

eig_active = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17, 18, 19, 20, 21, 22, 23];

t_final = 200;

N = length(t_plot); % Number of interpolated points
y_LPS = zeros(M, N);
y_LPSR = zeros(M, N);
y_IKK = zeros(M, N);
y_mRNA_R = zeros(M, N);
y_R = zeros(M, N);
y_NFkBn = zeros(M, N);
y_mRNA_IkBa = zeros(M, N);
y_IkBa = zeros(M, N);
y_P = zeros(M, N);
y_A = zeros(M, N);
y_E = zeros(M, N);
y_EPI = zeros(M, N);
y_cAMP = zeros(M, N);
y_M = zeros(M, N);
y_F = zeros(M, N);
y_mRNA_R_F = zeros(M, N);
y_R_F = zeros(M, N);
y_FR = zeros(M, N);
y_FRn = zeros(M, N);
y_R_EPI = zeros(M, N);
y_EPIR = zeros(M, N);
if eig
    e_d = zeros(M, N);
end

if init_parfor_progress
    parfor_progress(M);
end
parfor i=1:M
    s = struct();
    s.paramset = 'critical_transition';
    s.nonnegative = false;
    s.circadian = true;
    s.suppress_circadian = false;
    s.transition_logistic = true;
    s.t_pre = 408+24;
    s.t_final = t_final;
    s.eig = eig;
    % With ultradian disabled, 14, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
    % 35, and 36 aren't doing anything, which leaves the following:
    s.eig_active = eig_active;
    s.n_eig = 2000;
    s.t_ex = t_ex;

    LPS_0 = NaN;
    t_LPS = NaN;
    if (strcmp(vary, 't_LPS'))
        LPS_0 = LPS_0s;
        t_LPS = t_LPSs(i);
    elseif (strcmp(vary, 'LPS_0'))
        LPS_0 = LPS_0s(i);
        t_LPS = t_LPSs;
    end

    [t, y, e] = hem.model.run(LPS_0, t_LPS, s);

    if center_plots
        center_factor = t_LPS;
    else
        center_factor = floor(min(t_LPSs)/24)*24;
    end

    y_LPS(i,:) = interp1(t-s.t_pre-center_factor, y.LPS, t_plot);
    y_LPSR(i,:) = interp1(t-s.t_pre-center_factor, y.LPSR, t_plot);
    y_IKK(i,:) = interp1(t-s.t_pre-center_factor, y.IKK, t_plot);
    y_mRNA_R(i,:) = interp1(t-s.t_pre-center_factor, y.mRNA_R, t_plot);
    y_R(i,:) = interp1(t-s.t_pre-center_factor, y.R, t_plot);
    y_NFkBn(i,:) = interp1(t-s.t_pre-center_factor, y.NFkBn, t_plot);
    y_mRNA_IkBa(i,:) = interp1(t-s.t_pre-center_factor, y.mRNA_IkBa, t_plot);
    y_IkBa(i,:) = interp1(t-s.t_pre-center_factor, y.IkBa, t_plot);
    y_P(i,:) = interp1(t-s.t_pre-center_factor, y.P, t_plot);
    y_A(i,:) = interp1(t-s.t_pre-center_factor, y.A, t_plot);
    y_E(i,:) = interp1(t-s.t_pre-center_factor, y.E, t_plot);
    y_EPI(i,:) = interp1(t-s.t_pre-center_factor, y.EPI, t_plot);
    y_cAMP(i,:) = interp1(t-s.t_pre-center_factor, y.cAMP, t_plot);
    y_M(i,:) = interp1(t-s.t_pre-center_factor, y.M, t_plot);
    y_F(i,:) = interp1(t-s.t_pre-center_factor, y.F, t_plot);
    y_mRNA_R_F(i,:) = interp1(t-s.t_pre-center_factor, y.mRNA_R_F, t_plot);
    y_R_F(i,:) = interp1(t-s.t_pre-center_factor, y.R_F, t_plot);
    y_FR(i,:) = interp1(t-s.t_pre-center_factor, y.FR, t_plot);
    y_FRn(i,:) = interp1(t-s.t_pre-center_factor, y.FRn, t_plot);
    y_R_EPI(i,:) = interp1(t-s.t_pre-center_factor, y.R_EPI, t_plot);
    y_EPIR(i,:) = interp1(t-s.t_pre-center_factor, y.EPIR, t_plot);
    if eig
        e_d(i,:) = interp1(e.t-s.t_pre-center_factor, max(real(e.d), [], 2), t_plot);
    end
    parfor_progress;
end
if init_parfor_progress
    parfor_progress(0);
end

y_plot = struct();
for i=1:length(vars)
    y_plot.(vars{i}) = eval(sprintf('y_%s', vars{i}));
end
if eig
    y_plot.e_d = e_d;
end