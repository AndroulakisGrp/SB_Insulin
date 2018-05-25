function [t, y] = run(varargin)
global FF tt sigma G_star f_b K_f h_rho G_st tstart

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
p.addParamValue('glu', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
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
p.addParamValue('paramtest', 1, @(x)validateattributes(x, {'numeric'}, {'vector'}));

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
y0_insc=[0.1178 0.1 5 0.6 0.5 7.5115e-6 0 0 1e-5 1];
y0c=[y0_original' y0_met' y0_insc]';
y0_ins=[0.1178 0.1 5 0.6 0.5 7.5115e-6 0 0 1e-5 1 5];
y0=[y0_original' y0_met' y0_ins]';

% First, load all default parameter values
params_to_test = {'hh_circ', 'hh', 'hh_fen', 'hh_all', 'new', 'met', 'hh_jusko', 'hh_sig', 'hh_opt','glu'};
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


options = odeset('Nonnegative', 1:27, 'MaxStep',1.0e-2);

if ~s.stochastic
    %     keyboard;
    %tspan1  = [0 s.t_final];
    tspan1 = [0 1000];
   
    [t1, y1] = ode15s(@hem.model.model2, tspan1, y0 , options, s);
%    [t1, y1] = ode45(@hem.model.modelPan, tspan1, y0 , options, s);
    [tc, yc] = ode15s(@hem.model.modelc, tspan1, y0c , options, s);    
    t = t1;
    y = y1;
    
    
    assignin('base', 't1', t1)
    
    varsc = {'F', 'mRNA_R_F', 'R_F', 'FR' 'FRn', 'CRH', 'ACTH', 'F_per',...
        'percry', 'PERCRY', 'nucPERCRY', 'bmal', 'BMAL', 'nucBMAL',...
        'CLOCKBMAL','R_F_per', 'FR_per', 'FRn_per', 'M_F_per', 'FM_per',...
        'FMn_per', 'mRNA_P', 'P', 'mRNA_R_P','R_P', 'PR','P_cen',...
        'NAD','NAM','NMN','SIRT1','CLOCKBMALSIRT1','NAMPT','feed2','feed3','EntF',...
        'pgc1a','PGC1a','PGC1aN','PGC1aNa','FOXO1','Gluc','I','V','R','D','D_IR','Fins','gamma','rho','F2','G2'};
    vars = {'F', 'mRNA_R_F', 'R_F', 'FR' 'FRn', 'CRH', 'ACTH', 'F_per',...
        'percry', 'PERCRY', 'nucPERCRY', 'bmal', 'BMAL', 'nucBMAL',...
        'CLOCKBMAL','R_F_per', 'FR_per', 'FRn_per', 'M_F_per', 'FM_per',...
        'FMn_per', 'mRNA_P', 'P', 'mRNA_R_P','R_P', 'PR','P_cen',...
        'NAD','NAM','NMN','SIRT1','CLOCKBMALSIRT1','NAMPT','feed2','feed3','EntF',...
        'pgc1a','PGC1a','PGC1aN','PGC1aNa','FOXO1','Gluc','I','V','R','D','D_IR','Fins','gamma','rho','F2','G2','b_V'};
    ids=1:length(vars);
   
    y_old = y;
    y_oldc = yc;
    y = struct;
    yc = struct;
    for j=1:length(ids)
        id = ids(j);
        y.(vars{j}) = y_old(:,id)';
    end
    for j=1:(length(ids)-1)
        id = ids(j);
        yc.(varsc{j}) = y_oldc(:,id)';
    end
    
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
        'NAD','NAM','NMN','SIRT1','CLOCKBMALSIRT1','NAMPT','feed2','feed3','EntF',...
        'pgc1a','PGC1a','PGC1aN','PGC1aNa','FOXO1','Gluc','I','V','R','D','D_IR','Fins','gamma','rho','F2','G2','b_V'};
    ids=1:length(vars);
    y_old = y;
    y = struct;
    for j=1:length(ids)
        id = ids(j);
        y.(vars{j}) = y_old(:,id)';
    end
end
% keyboard;



%% Plotting
if s.make_plots
 %   vars_to_plot = {'F','PERCRY','CLOCKBMAL','PGC1aN','Gluc'};
%    titles = {'Cortisol','PER/CRY','CLOCK/BMAL1','PGC1-\alpha(N)','Pck1/G6pc'};
    vars_to_plot = {'F','BMAL','PERCRY','SIRT1','FOXO1','Gluc'};
    titles = {'Cortisol','BMAL1','PER/CRY','SIRT1','FOXO1','Pck1/G6pc'};
    
    n = length(vars_to_plot);
    
    %     Make appropriately sized subplot
    nx = ceil(sqrt(n));
    ny = ceil(sqrt(n));
    if nx*(ny-1) >= n
        ny = ny - 1;
    end
    

    for i=1:n
        var = vars_to_plot{i};
        hold on
        subplot(2,3,i);
        hold on

        plot(t, mean(y.(var), 1),'k:','LineWidth',2);

        
        xlim([480 480+24]);
        set(gca,'XTick',480:12:480+24);
        set(gca,'XTickLabel',{'12am','12pm','12am'});
        hold on;
        title(titles(i),'FontSize',12);
    end
    
    time=linspace(480,480+(24*21),3000);
    lightcontrol=zeros(length(time),1);   lighttest=zeros(length(time),1);
    for i=1:length(time)
        if (mod(time(i)+s.lst3,24.0000) > s.sp)
            lightcontrol(i)=1;
        else
            lightcontrol(i)=0;
        end
    end

end

if s.extra_plots
    % Plot all variables along with their max and min values
    %         figure;
    vars = {'F', 'mRNA_R_F', 'R_F', 'FR' 'FRn', 'CRH', 'ACTH', 'F_per',...
        'percry', 'PERCRY', 'nucPERCRY', 'bmal', 'BMAL', 'nucBMAL',...
        'CLOCKBMAL','R_F_per', 'FR_per', 'FRn_per', 'M_F_per', 'FM_per',...
        'FMn_per', 'mRNA_P', 'P', 'mRNA_R_P','R_P', 'PR','P_cen',...
        'NAD','NAM','NMN','SIRT1','CLOCKBMALSIRT1','NAMPT','feed2','feed3','EntF',...
        'pgc1a','PGC1a','PGC1aN','PGC1aNa','FOXO1','Gluc'};
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
        plot(t, mean(y.(var), 1),'b','LineWidth',2);
        hold on
        plot(tc, mean(yc.(var),1),'k:','LineWidth',2);
        xlim([840 840+24]);
        set(gca, 'XTick', 840:12:840+24);
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

    I_0=9.2586397e-9;   % ug
    N_c=1000;    N=2.76e6;  % # of pancreatic beta-cells
    Finsc=mean(yc.Fins,1);      Fins=mean(y.Fins,1);
    G=zeros(1,length(t));       f=zeros(1,length(t));
    Gc=zeros(1,length(tc));       fc=zeros(1,length(tc));
    toc
    
    for i=1:length(tc)   
        if tc(i) < 6
            Gc(i)=1;
        elseif mod(tc(i)-6,24) <= 0
            Gc(i)=1;
        elseif mod(tc(i)-6,24) < 1
            Gc(i)=G_st;
        elseif mod(tc(i)-6,24) <= 5
            Gc(i)=1;
        elseif mod(tc(i)-6,24) < 6
            Gc(i)=G_st;
        elseif mod(tc(i)-6,24) <= 10
            Gc(i)=1;
        elseif mod(tc(i)-6,24) < 11
            Gc(i)=G_st;
        elseif mod(tc(i)-6,24) >= 11
            Gc(i)=1;
        end
    end
    
    for i=1:length(t)   
        if t(i) < tstart
            G(i)=1;
        elseif mod(t(i)-tstart,24) <= 0
            G(i)=1;
        elseif mod(t(i)-tstart,24) < 1
            G(i)=G_st*0.5;
        elseif mod(t(i)-tstart,24) <= 5
            G(i)=1;
        elseif mod(t(i)-tstart,24) < 6
            G(i)=G_st*0.75;
        elseif mod(t(i)-tstart,24) <= 10
            G(i)=1;
        elseif mod(t(i)-tstart,24) < 11
            G(i)=G_st;
        elseif mod(t(i)-tstart,24) >= 11
            G(i)=1;
        end
    end
    
    for i=1:length(t)
        if G(i) < G_star
            f(i)=f_b;
        else
            f(i)=f_b+(1-f_b)*((G(i)-G_star)/(K_f+G(i)-G_star));
        end
    end
    
    %Z = rho_bar ./ a .*(b*T + (D_IRb - b./a) .* (1-exp(-a.*T)));
    %IS = I_0 .* Z .* f .*N;
    ISRc = I_0 .* sigma .* Finsc .*  fc .* N;     % ug/hr
    ISR = I_0 .* sigma .* Fins .*  f .* N;
    
    figure
    vars={'I','V','R','D','D_IR','Fins','gamma','rho','F2','G2'};
    varname={'I','V','R','D','D_{IR}','Fins','\gamma','\rho','F2','G2'};
    for i=1:length(vars)
        var=vars{i};
        hold on
        subplot(4,4,i)
        hold on
        plot(tc,mean(yc.(var), 1))
        plot(t,mean(y.(var),1))
        xlim([960,960+24])
        xlabel('time(h)'); title(varname(i))
        hold on
    end
    
    subplot(4,4,length(varname)+1)
    hold on
    plot(tc,ISRc)
    plot(t,ISR)
    xlim([960,960+24])
    xlabel('time(h)');    title('ISR');
    
    subplot(4,4,length(varname)+2)
    hold on
    plot(tc,Gc)
    plot(t,G)
    xlim([960,960+24])
    xlabel('time(h)');    title('G');
    
    subplot(4,4,length(varname)+3)
    hold on
    plot(tc, fc)
    plot(t,f)
    xlim([960,960+24])
    xlabel('time(h)');    title('f');
    hold on
    
    name=num2str(tstart);
    name2=strcat('bigdinnerISR_', name);
    save(name2, 'ISR')

    
end
%%
if s.ss_plots
     vars = {'F', 'mRNA_R_F', 'R_F', 'FR' 'FRn', 'CRH', 'ACTH', 'F_per',...
         'percry', 'PERCRY', 'nucPERCRY', 'bmal', 'BMAL', 'nucBMAL',...
         'CLOCKBMAL','R_F_per', 'FR_per', 'FRn_per', 'M_F_per', 'FM_per',...
         'FMn_per', 'mRNA_P', 'P', 'mRNA_R_P','R_P', 'PR','P_cen',...
         'NAD','NAM','NMN','SIRT1','CLOCKBMALSIRT1','NAMPT','feed2','feed3','EntF',...
        'pgc1a','PGC1a','PGC1aN','PGC1aNa','FOXO1','Gluc'};
    %vars = {'F','PERCRY'};
    %titles = {'Cortisol','PERCRY'};

n = length(vars);
    figure
    hold on
    for i=1:n
        var = vars{i};
        subplot(ceil(sqrt(n)), ceil(sqrt(n)), i);
        plot(t, mean(y.(var), 1),'k-');
  %      hold on
 %       plot(tc, mean(yc.(var),1),'k:');
%         hold on
%         plot(t, mean(y.(var), 1),'b--');
        %xlim([720 720+240]);
        %xlim([500-240 2000])
        xlim([400 700]);
        %set(gca, 'XTick', 720:24:720+240);
        %set(gca, 'XTickLabel',{'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'})
        title(vars(i),'FontSize',10);
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
