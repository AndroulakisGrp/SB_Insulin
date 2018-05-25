function [t, y] = run(varargin)
%HEM.MODEL.RUN Run a simulation of the human endotoxemia model.
%   [T,Y] = HEM.MODEL.RUN(LPS,T_LPS) runs a simulation for a dose of
%   endotoxin (LPS) given at a certain time (t_LPS). Y is a struct
%   containing all of the model variables, and T is the corresponding timeare
%   vector.
%
%   The arguments LPS,T_LPS can be followed by parameter/value paris to
%   specify additional options (these can also be passed in a struct):
%
%      OUTPUT/ANALYSIS
%      'eig' is a boolean defining whether or not the eigenvalues of the
%         Jacobian matrix should be estimated. If true, then additional
%         output E is returned, giving the time (E.t), eigenvalues (E.d),
%         and eigenvectors (E.v). At times prior to t_pre (see below), E.d
%         and E.v are 0s.
%      'eig_active' is a list of variable indices (as in HEM.MODEL.MODEL)
%         that will be used to calculate the Jacobian. This allows you to
%         ignore irrelevant variables, e.g. ignoring ultradian variables
%         if you have ultradian rhythms disabled.
%      'n_eig' approximate number of eigenvalues to calculate, evenly
%         spaced through time. If this is too low, eigenvalues and
%         eigenvectors may be poorly traced through time.
%
%      EXOGENOUS INFUSION OR CHRONIC PRODUCTION
%      'F_ex' defines exogenous cortisol infusion into the dF/dt equation
%      'EPI_ex' defines exogenous epinephrine infusion into the dEPI/dt
%         equation
%      'LPS_ex' defines exogenous LPS infusion into the dLPS/dt equation
%      't_ex' is a vector containing the start and stop times of infusion
%      'F_chronic' defines exogenous cortisol infusion into the dF_ex/dt
%         equation
%      'P_chronic' defines chronic production of P
%
%      INTEGRATION
%      'nonnegative' is a boolean for whether to enable the NonNegative
%         option for ode45 or not (default true). IF YOU SET THIS TO FALSE,
%         USE THE NONNEGATIVE PARAMSET or some other paramset that accoutns
%         for this.
%      't_final' is the final time point in the integration
%      't_pre' is how long to run the integration (letting it reach steady
%         state) before doing anything else; this only matters for the
%         circadian model, generally
%      'stochastic_ensemble' is a boolean; true: stochastic integration
%         with the Euler-Maruyama method; false: deterministic integration
%         using ODE45 (default is false/deterministic); THE FOLLOWING
%         INTEGRATION SETTINGS ONLY APPLY IF stochastic_ensemble IS true
%      'noise_type' can be 'additive', 'multiplicative', or 'langevin';
%         anything besides 'additive' (the default) is not well tested
%      'h' is the time step for stochastic integration; the default value
%         of 0.0234 should generally work well unless you make some large
%         changes to the model
%      'mu' is the noise level; for additive noise, it is multiplied by the
%         average value of the variable
%      'n_cells' is the number of cells (leukocytes) to simulate; anything
%         besides 1 is not well tested
%
%      MODEL
%      'circadian' is a boolean (default true) that can enable/disable the
%         circadian model components
%      'suppress_circadian' is a boolean (default false) that can
%         enable/disable the suppression of circadian rhythms at high
%         levels of inflammatory mediators.
%      'pulsatile' is a boolean (default false) that can enable/disable the
%         pulsatile cortisol model
%      'NFkB_spiky' is a boolean (default false) that can enable an
%         alternative model of NFkB; not well tested
%      'com_mult' is a multiplier on some parameters linking model
%         subcompartments;  not well defined, don't use it without checking
%         what it does
%      'P_inhibits_FRn' is a boolean indicating if P should inhibit
%         translocation of FR to the nucleus (not well tested)
%      'F_noise_IC_factor' is a factor that is multiplied by the noise
%         amplitude in certain variables; not well defined, don't use it
%         without checking what it does
%      'transition_logistic' uses the logistic growth equation for LPS, as
%         defined in the transition IEEE paper. With this enabled, the
%         second value of t_ex defines the time at which LPS elimination
%         begins by increasing the degradation rate.
%      'paramset' is a string corresponding to a subfolder of +hem/params;
%         is set, parameters in this folder replace default parameter
%         values
%      Finally, each parameter file (hh, hh_all, hh_circ, hh_fen, hh_jusko,
%         hh_lps, hh_sig) will be replaced by a vector of the same size
%         if one is passed here with that name
%
%      PLOTTING
%      'make_plots' is a boolean that enables/disables the ploting of
%          pre-defined variables.
%      'plot_data' is a boolean that enables/disables the plot of the
%          experimental data that has been used in the JTB paper
%      'extra_plots' is a boolean that enables/disables the plot of
%          additional plots when make_plots is true that are also
%          pre-defined.
%
%   See also HEM.PLOT.VARS, HEM.MODEL.MODEL, HEM.MODEL.IPFM.

%% Function parameters
p = inputParser;
p.StructExpand = true;


% INTEGRATION
% Enforce nonnegativity?
p.addParamValue('nonnegative', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% How long to integrate?
p.addParamValue('t_final', 48, @(x)validateattributes(x, {'numeric'}, {'scalar'}));
% How long to run before starting everything? This is helpful if parameters
% or equations are changing, so initial conditions might also change.
p.addParamValue('t_pre', 48, @(x)validateattributes(x, {'numeric'}, {'scalar'}));
% false: deterministic, one cell; true: stochastic, multiple cells
% Time range of exogenous LPS dosing
p.addParamValue('t_ex_LPS', [0, 0], @(x)validateattributes(x, {'numeric'}, {'vector', 'size', [1, 2]}));
% Stochastic?
p.addParamValue('stochastic', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));


% MODEL
% Include circadian components?
p.addParamValue('circadian', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% Use LPS logistic growth term?
p.addParamValue('transition_logistic', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% Use HPA model of cortisol?
p.addParamValue('HPA', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% Load parameters from files if they are not passed explicitly
p.addParamValue('hh_circ', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('hh', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('hh_fen', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('hh_all', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('new', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('hh_jusko', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('hh_sig', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('hh_opt', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('met', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
% Period of scotoperiod
p.addParamValue('sp', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
% Light stimulus minimum and maximum
p.addParamValue('lst1', 1, @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('lst2', 0, @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('lst3', 6, @(x)validateattributes(x, {'numeric'}, {'vector'}));


% PLOTTING
% Show any plots?
p.addParamValue('make_plots', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% Show data (for baseline conditions) on main plot?
p.addParamValue('plot_data', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% If make_plots is also true, show some additional plots (currently 2: all
% of the variables in subplots and individual vs average values of NFkBn)
p.addParamValue('extra_plots', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
p.addParamValue('ss_plots', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% Color for the main plot
p.addParamValue('color', 'k', @(x)validateattributes(x, {'char'}, {}));
% Line style for the main plot
p.addParamValue('line_style', '-', @(x)validateattributes(x, {'char'}, {}));

p.addParamValue('runx', 0, @(x)validateattributes(x, {'numeric'}, {'scalar'}));
p.addParamValue('runy', 0, @(x)validateattributes(x, {'numeric'}, {'scalar'}));
p.addParamValue('fdpd', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));

% Food stimulus minimum and maximum
p.addParamValue('fdhigh', 1, @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('fdlow', 0, @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('fdphase', 6, @(x)validateattributes(x, {'numeric'}, {'vector'}));

p.addParamValue('cells', 1, @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('count', 1, @(x)validateattributes(x, {'numeric'}, {'vector'}));

p.parse(varargin{:});
s = p.Results;

s.t_final = s.t_final + s.t_pre;

%% Parameter and initial condition defaults

f = hem.util.get_hem_folder; % Path to +hem folder

% Initial conditions - should be recalculated if parameters or equations
% are changed
% y0 = [s.LPS, 0, 0, 0.93222, 0.9163, 1.3913e-05, 1.9434e-10, 0.1233, 0.96769, 1.2401, 0.86885, 0.99495, 1, 0.95934, 0.46334, 0, 0.98761, 25.251, 522.98, 0, 0.016839, 1.0013, 0, 0, 1, 1, 1, 1, 0, 25.251, 522.98, 0, 0.016839, 0];
y0_original = importdata(sprintf('%s/param/ICs.txt', f));
y0_met = importdata(sprintf('%s/param/ICmet.txt', f));
y0=[y0_original' y0_met']';
% f = hem.util.get_hem_folder;
% incon = importdata(sprintf('%s/param/hh_opt.txt', f));
% y0(1) = random('Uniform',0,1);
% y0(2) = random('Uniform',0,1);
% y0(3) = random('Uniform',0.5,1);
% y0(4)=  random('Uniform',0.5,1);
% y0(5)=  random('Uniform',0.5,1);
% y0(6)=  random('Uniform',0.5,1);
% y0(7)=  random('Uniform',0.5,1);
% y0(8)=  random('Uniform',0.5,1);
% y0(9)= random('Uniform',0.5,1);
% y0(10)=  random('Uniform',0.5,1);
% y0(11)=  random('Uniform',0.5,1);
% y0(12)= random('Uniform',0.5,1);
% y0(13)=  random('Uniform',0.5,1);
% y0(14)=  random('Uniform',0.5,1);
% y0(15)= random('Uniform',0.5,1);
% y0(16)= random('Uniform',0.5,1);
% y0(17)= random('Uniform',0.5,1);
% y0(18)= random('Uniform',0.5,1);
% y0(19)= random('Uniform',0.5,1);
% y0(20)= random('Uniform',0.5,1);
% y0(21)= random('Uniform',0.5,1);
% y0(22)= incon(80);
% y0(23)= incon(81);
% y0(24)= incon(82);
% y0(25)= incon(83);
% y0(26)= incon(84);
% y0(27)= incon(85);
% y0(28)= incon(86);
% y0(29)= incon(87);
% y0(30)= incon(88);


% First, load all default parameter values
params_to_test = {'hh_circ', 'hh', 'hh_fen', 'hh_all', 'new', 'met', 'hh_jusko', 'hh_sig', 'hh_opt'};
load_param = zeros(size(params_to_test)); % Was parameter passed to function or loaded from file?
for i=1:length(params_to_test)
    if isempty(s.(params_to_test{i}))
        s.(params_to_test{i}) = importdata(sprintf('%s/param/%s.txt', f, params_to_test{i}));
        load_param(i) = true;
    else
        load_param(i) = false;
    end
end

% Is circadian enabled or not?
if ~s.circadian
    s.hh_circ = zeros(size(s.hh_circ));
end


%% Integration


options = odeset('Nonnegative', 1:27);

if ~s.stochastic
    %     keyboard;
    tspan1  = [0 2016];
    
    [t1, y1] = ode45(@hem.model.cortisol_only, tspan1, y0 , options, s);
    [tc, yc] = ode45(@hem.model.modelc, tspan1, y0 , options, s);    
    t = t1;
    y = y1;
    
    assignin('base', 't1', t1)
    
    
    vars = {'F', 'mRNA_R_F', 'R_F', 'FR' 'FRn', 'CRH', 'ACTH', 'F_per',...
        'percry', 'PERCRY', 'nucPERCRY', 'bmal', 'BMAL', 'nucBMAL',...
        'CLOCKBMAL','R_F_per', 'FR_per', 'FRn_per', 'M_F_per', 'FM_per',...
        'FMn_per', 'mRNA_P', 'P', 'mRNA_R_P','R_P', 'PR','P_cen',...
        'NAD','NAM','NMN','SIRT1','CLOCKBMALSIRT1','NAMPT','feed2','feed3','EntF'};
    ids=1:length(vars);
   
    y_old = y;
    y_oldc = yc;
    y = struct;
    yc = struct;
    for j=1:length(ids)
        id = ids(j);
        y.(vars{j}) = y_old(:,id)';
        yc.(vars{j}) = y_oldc(:,id)';
    end
    
    % for tim=1:length(t)
    % % %     if t(tim)>s.t_pre+s.lst && t(tim)<s.t_pre+s.lst+1
    % % %         light(tim)=1;
    % % Simulate regular step light
    % if (mod(t(tim)+s.lst3,24) >s.sp)
    %     light(tim)=s.lst1;
    % else
    %     light(tim)=s.lst2;
    % end
    %
    % % % Simulate increasing light (maximum at the center)
    % %     if (mod(t(tim)+6,24.0000) >=12 && mod(t(tim)+6,24.0000) <14)
    % %         light(tim)=0.4;
    % %     else if (mod(t(tim)+6,24.0000) >=14 && mod(t(tim)+6,24.0000) <16)
    % %         light(tim)=0.6;
    % %         else if (mod(t(tim)+6,24.0000) >= 16 && mod(t(tim)+6,24.0000) < 20)
    % %         light(tim)=1;
    % %             else if(mod(t(tim)+6,24.0000) >= 20 && mod(t(tim)+6,24.0000) < 22)
    % %         light(tim)=0.6;
    % %                 else if (mod(t(tim)+6,24.0000) >= 22 && mod(t(tim)+6,24.0000) <= 24)
    % %         light(tim)=0.4;
    % %                     else if (mod(t(tim)+6,24.0000) <12)
    % %         light(tim)=0;
    % %                     end
    % %                 end
    % %             end
    % %         end
    % %     end
    % % end
    %
    % % % Simulate increasing light (maximum at the left side)
    % %     if (mod(t(tim)+6,24.0000) >=12 && mod(t(tim)+6,24.0000) <14)
    % %         light(tim)=0.4;
    % %     else if (mod(t(tim)+6,24.0000) >=14 && mod(t(tim)+6,24.0000) <18)
    % %         light(tim)=1;
    % %         else if (mod(t(tim)+6,24.0000) >= 18 && mod(t(tim)+6,24.0000) < 20)
    % %         light(tim)=0.6;
    % %             else if(mod(t(tim)+6,24.0000) >= 20 && mod(t(tim)+6,24.0000) < 22)
    % %         light(tim)=0.6;
    % %                 else if (mod(t(tim)+6,24.0000) >= 22 && mod(t(tim)+6,24.0000) <= 24)
    % %         light(tim)=0.4;
    % %                     else if (mod(t(tim)+6,24.0000) <12)
    % %         light(tim)=0;
    % %                     end
    % %                 end
    % %             end
    % %         end
    % %     end
    % % end
    %
    %
    % % % Simulate increasing light (maximum at the right side)
    % %     if (mod(t(tim)+6,24.0000) >=12 && mod(t(tim)+6,24.0000) <14)
    % %         light(tim)=0.4;
    % %     else if (mod(t(tim)+6,24.0000) >=14 && mod(t(tim)+6,24.0000) <16)
    % %         light(tim)=0.6;
    % %         else if (mod(t(tim)+6,24.0000) >= 16 && mod(t(tim)+6,24.0000) < 18)
    % %         light(tim)=0.6;
    % %             else if(mod(t(tim)+6,24.0000) >= 18 && mod(t(tim)+6,24.0000) < 22)
    % %         light(tim)=1;
    % %                 else if (mod(t(tim)+6,24.0000) >= 22 && mod(t(tim)+6,24.0000) <= 24)
    % %         light(tim)=0.4;
    % %                     else if (mod(t(tim)+6,24.0000) <12)
    % %         light(tim)=0;
    % %                     end
    % %                 end
    % %             end
    % %         end
    % %     end
    % % end
    %
    % % % Simulate jet lag at time t=2400hr
    % % if t(tim)>2400
    % %     if (mod(t(tim)+s.lst3,24) >s.sp)
    % %         light(tim)=s.lst1;
    % %     else
    % %         light(tim)=s.lst2;
    % %     end
    % % else
    % %     if (mod(t(tim)+s.lst3,24) >s.sp)
    % %         light(tim)=s.lst1;
    % %     else
    % %         light(tim)=s.lst2;
    % %     end
    % % end
    %
    % % % Experiment of Jun et al. to play
    % %         if t(tim)>s.t_pre+19 && t(tim)<s.t_pre+19+6.7
    % %             light(tim)=67;
    % %             % else if t>=s.t_pre+24 && t<s.t_pre+96+7*24+19+6.7
    % %             %         light=0.02;
    % %             %     end
    % %         end
    %
    %
    % end
    % keyboard
    
else
    % Stochastic constant noise integration
    y0=y0';
    mu = 0.03*ones(size(y0));
    M = length(y0);                    % num of variables
    dt = 0.1; N = s.t_final/dt;
    dW = sqrt(dt)*randn(N, M);         % Brownian increments
    y_temp = y0;
    y = [];
    y(1,:) = y_temp;
    %     keyboard;
    for j = 1:N
        t = dt*j;
        dy = hem.model.model(t, y_temp,s);
        y_temp = y_temp + dt*dy + mu.*dW(j,:);
        y(j+1,:) = y_temp;
        %        keyboard;
    end
    t = 0:dt:s.t_final;
    
    % keyboard;
    vars = {'F', 'mRNA_R_F', 'R_F', 'FR' 'FRn', 'CRH', 'ACTH', 'F_per',...
        'percry', 'PERCRY', 'nucPERCRY', 'bmal', 'BMAL', 'nucBMAL',...
        'CLOCKBMAL','R_F_per', 'FR_per', 'FRn_per', 'M_F_per', 'FM_per',...
        'FMn_per', 'mRNA_P', 'P', 'mRNA_R_P','R_P', 'PR','P_cen',...
        'NAD','NAM','NMN','SIRT1','CLOCKBMALSIRT1','NAMPT','feed2','feed3','EntF'};
    ids=1:length(vars);
    y_old = y;
    y = struct;
    for j=1:length(ids)
        id = ids(j);
        y.(vars{j}) = y_old(:,id)';
    end
end
% keyboard;

% Set light and Feeding
light = zeros(1,length(t1));
for i=1:length(t1)
    if (mod(t1(i)+s.lst3,24.0000) > s.sp)
        light(i)=s.lst1;
    else
        light(i)=s.lst2;
    end
end

feed = zeros(1,length(t1));
for i=1:length(t1)
   if t1(i)<720
        if (mod(t1(i)-s.fdphase,24) < s.fdpd)
            feed(i)=s.fdhigh;
        else
            feed(i)=s.fdlow;
        end
    else
        if (mod(t1(i)-s.fdphase,24) < s.fdpd)      % reverse feeding
            feed(i)=s.fdlow;
        else
            feed(i)=s.fdhigh;
        end
    end
end

%% Plotting
if s.make_plots
    
    vars_to_plot = {'CRH','F_per'};    %'CRH','ACTH','F','P','A','EPI'
    titles = {'CRH','F_{periphery}'};
    
    n = length(vars_to_plot);
    
    %     Make appropriately sized subplot
    nx = ceil(sqrt(n));
    ny = ceil(sqrt(n));
    if nx*(ny-1) >= n
        ny = ny - 1;
    end
    
    for i=1:n
        
        subplot(2, 2, i+2);
        
        if s.plot_data
            [data_t, data_y] = hem.model.get_data(vars_to_plot{i});
            data_t = [data_t, 24+data_t, 48+data_t]; %#ok<AGROW>
            data_y = [data_y, data_y, data_y]; %#ok<AGROW>
        else
            data_t = [];
            data_y = [];
        end
        
        if strcmp('Fboth', vars_to_plot{i})
            y_plot = y.F + y.F_ex;
        else
            y_plot = y.(vars_to_plot{i});
        end
        %         rectangle('Position',[18,0,12,2],'FaceColor',[0.8 0.8 0.8])
        hem.model.plot_variable(t, y_plot, data_t, data_y, s.t_pre,s.line_style,s.color);
        ylabel(titles{i});
        
        set(gca, 'XTick', [0:12:24]);
        set(gca, 'XTickLabel',{'12am', '12pm', '12am'}) % {'1st week', '2nd week','3rd week', '4th week'}   {'12am', '12am'}
        xlim([0 24]);
        grid on;
        %         rectangle('Position',[18,0,12,1.5],'FaceColor',[0.5 0.5 0.5])
        %         suplabel('Time (hr)', 'x', [.1 .1 .85 .85]);
        
    vars_to_plot = {'F','percry','bmal','CLOCKBMAL','SIRT1','NAD','NAMPT'};
    titles = {'Cortisol','PerCry_{mRNA}','Bmal_{mRNA}','CLOCK/BMAL1',...
        'SIRT1','NAD^+','NAMPT'};
    
    n = length(vars_to_plot);
    
    %     Make appropriately sized subplot
    nx = ceil(sqrt(n));
    ny = ceil(sqrt(n));
    if nx*(ny-1) >= n
        ny = ny - 1;
    end
    
    figure
    hold on
    
    for i=1:n
        
        subplot(3, 3, i+2);
        
        if s.plot_data
            [data_t, data_y] = hem.model.get_data(vars_to_plot{i});
            data_t = [data_t, 24+data_t, 48+data_t]; %#ok<AGROW>
            data_y = [data_y, data_y, data_y]; %#ok<AGROW>
        else
            data_t = [];
            data_y = [];
        end
        
        if strcmp('Fboth', vars_to_plot{i})
            y_plot = y.F + y.F_ex;
        else
            y_plot = y.(vars_to_plot{i});
        end
        %         rectangle('Position',[18,0,12,2],'FaceColor',[0.8 0.8 0.8])
        hem.model.plot_variable(t, y_plot, data_t, data_y, s.t_pre,s.line_style,s.color);
        ylabel(titles{i});
        
        set(gca, 'XTick', [0:12:24]);
        set(gca, 'XTickLabel',{'12am', '12pm', '12am'}) % {'1st week', '2nd week','3rd week', '4th week'}   {'12am', '12am'}
        xlim([0 24]);
        grid on;
    end
    
    subplot(3,3,1)
    hold on
    plot(t,light,s.color)
    xlim([s.t_pre s.t_pre+24]);
    ylim([0 1.2])
    grid on
    ylabel('Light');
    set(gca, 'XTick', s.t_pre:12:s.t_pre+24);
    set(gca, 'XTickLabel',{'12am', '12pm', '12am'})
    
    subplot(3,3,2)
    hold on
    plot(t,feed,s.color)
    xlim([s.t_pre s.t_pre+24]);
    ylim([0 1.2])
    grid on
    ylabel('Feeding');
    set(gca, 'XTick', s.t_pre:12:s.t_pre+24);
    set(gca, 'XTickLabel',{'12am', '12pm', '12am'})
    %     hold on;
    %     subplot(2,2,1)
    %     hold on
    %     plot(t,light,'Color',s.color)
    %     xlim([s.t_pre s.t_pre+48]);
    %     grid on
    %     ylabel('Light')
    %     set(gca, 'XTick', [s.t_pre:24:s.t_pre+48]);
    %     set(gca, 'XTickLabel',{'12am', '12am'}) % {'1st week', '2nd week','3rd week', '4th week'}   {'12am', '12am'}
    
    
    %     % Calculate the period of light
    %     s.t_pre=480;
    %     s.t_final=800;
    %     Fs = 10;    % Sampling frequency that will be used in the FFT.
    %     b = (0:1/Fs:(s.t_final))+s.t_pre;   % Interpolation of signal in order to be used in FFT.
    %     e3 = interp1(t,light,b);
    %     e3 = e3 - mean(e3);
    %     NFFT = 20000;
    %     E3 = fft(e3,NFFT);
    %     m=max(E3);
    %     ang = angle(m);
    %     f = Fs*linspace(0,1,NFFT);
    %     [o3,l3] = max(E3(1:NFFT/2));
    %     T3 = f(l3)
end

if s.extra_plots
    % Plot all variables along with their max and min values
    %         figure;
    vars = {'F', 'mRNA_R_F', 'R_F', 'FR' 'FRn', 'CRH', 'ACTH', 'F_per',...
        'percry', 'PERCRY', 'nucPERCRY', 'bmal', 'BMAL', 'nucBMAL',...
        'CLOCKBMAL','R_F_per', 'FR_per', 'FRn_per', 'M_F_per', 'FM_per',...
        'FMn_per', 'mRNA_P', 'P', 'mRNA_R_P','R_P', 'PR','P_cen',...
        'NAD','NAM','NMN','SIRT1','CLOCKBMALSIRT1','NAMPT','feed2','feed3','EntF'};
    n = length(vars);
%     figure
%     hold on
%     for i=1:n
%         var = vars{i};
%         subplot(ceil(sqrt(n)), ceil(sqrt(n)), i);
%         plot(t, mean(y.(var), 1),s.color);
%         xlim([s.t_pre s.t_pre+24]);
%         set(gca, 'XTick', s.t_pre:12:s.t_pre+24);
%         set(gca, 'XTickLabel',{'12am', '12pm', '12am'})
%         hold on;
%         title(var,'FontSize',10);
%     end
%     subplot(ceil(sqrt(n)), ceil(sqrt(n)), i+1);
%     plot(t, feed,s.color);
%     xlim([s.t_pre s.t_pre+24]);
%     set(gca, 'XTick', s.t_pre:12:s.t_pre+24);
%     set(gca, 'XTickLabel',{'12am', '12pm', '12am'})
%     hold on;
%     title('feed','FontSize',10);
%     ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off',...
%         'Visible','off','Units','normalized', 'clipping' , 'off');
%     text(0.5, 1,'\bf Feeding from 6am-6pm','HorizontalAlignment',...
%         'center','VerticalAlignment', 'top')
    
    figure
    hold on
    for i=1:n
        var = vars{i};
        subplot(ceil(sqrt(n)), ceil(sqrt(n)), i);
        plot(t, mean(y.(var), 1),s.color);
        hold on
        plot(tc, mean(yc.(var),1),'k:');
        xlim([s.t_pre+960 s.t_pre+24+960]);
        set(gca, 'XTick', s.t_pre+960:12:s.t_pre+24+960);
        set(gca, 'XTickLabel',{'12am', '12pm', '12am'})
        hold on;
        title(var,'FontSize',10);
    end
%     subplot(ceil(sqrt(n)), ceil(sqrt(n)), i+1);
%     plot(t, feed,s.color);
%     xlim([s.t_pre+480 s.t_pre+24+480]);
%     set(gca, 'XTick', s.t_pre+480:12:s.t_pre+24+480);
%     set(gca, 'XTickLabel',{'12am', '12pm', '12am'})
%     hold on;
%     title('feed','FontSize',10);
%     ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off',...
%         'Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,'\bf Reversed Feeding','HorizontalAlignment',...
        'center','VerticalAlignment', 'top')

end
%%
if s.ss_plots
feed = zeros(1,length(t1));
for i=1:length(t1)
   if t1(i)<720
        if (mod(t1(i)-s.fdphase,24) < s.fdpd)
            feed(i)=s.fdhigh;
        else
            feed(i)=s.fdlow;
        end
    else
        if (mod(t1(i)-s.fdphase,24) < s.fdpd)      % reverse feeding
            feed(i)=s.fdlow;
        else
            feed(i)=s.fdhigh;
        end
    end
end
% i=1;
% feed=zeros(1,length(t));
% for i=1:length(t)
% if (mod(t(i)-s.fdphase,24) < s.fdpd)
%     feed(i)=s.fdhigh;
% else
%     feed(i)=s.fdlow;
% end
% end
i=1;
light = zeros(1,length(t1));
for i=1:length(t1)
    if (mod(t1(i)+s.lst3,24.0000) > s.sp)
        light(i)=s.lst1;
    else
        light(i)=s.lst2;
    end
end
%     vars = {'F', 'mRNA_R_F', 'R_F', 'FR' 'FRn', 'CRH', 'ACTH', 'F_per',...
%         'percry', 'PERCRY', 'nucPERCRY', 'bmal', 'BMAL', 'nucBMAL',...
%         'CLOCKBMAL','R_F_per', 'FR_per', 'FRn_per', 'M_F_per', 'FM_per',...
%         'FMn_per', 'mRNA_P', 'P', 'mRNA_R_P','R_P', 'PR','P_cen',...
%         'NAD','NAM','NMN','SIRT1','CLOCKBMALSIRT1','NAMPT','feed2','feed3','EntF'};
    vars = {'F','CRH','PERCRY','BMAL','NAD','SIRT1'};

n = length(vars);
    figure
    hold on
    i=1;
    for i=1:n
        var = vars{i};
        subplot(ceil(sqrt(n)), ceil(sqrt(n)), i);
        plot(t, mean(y.(var), 1),'k-');
        hold on
        plot(tc, mean(yc.(var),1),'k:');
%         hold on
%         plot(t, mean(y.(var), 1),'b--');
xlim([0 2000]);
        %xlim([720-240 720+480]);
        %xlim([140 300]);
%         set(gca, 'XTick', 960:12:984);
%         set(gca, 'XTickLabel',{'12am', '12pm', '12am'})
        title(var,'FontSize',10);
    end
%     hold on
%     subplot(3,3,1);
%     plot(t1, light ,'k-');
%     xlim([960 984]);
%     ylim([0 1.2]);
%     set(gca, 'XTick', 960:12:984);
%     set(gca, 'XTickLabel',{'12am', '12pm', '12am'})
%     title('Light','FontSize',10);
%     hold on
%     subplot(3,3,2);
%     plot(t, feed ,'b--');
%     hold on
%     plot(tr, feedr, 'k-');
%     xlim([960 984]);
%     ylim([0 1.2]);  
%     set(gca, 'XTick', 960:12:984);
%     set(gca, 'XTickLabel',{'12am', '12pm', '12am'})
%     title('Feed','FontSize',10);
end
end
