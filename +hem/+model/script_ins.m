close all;
clear;
clc;


tic
figure(1)
figure(2) 

for k=1:12
    
    
    global tt FF sigma G_star f_b K_f h_rho G_st tstart
    G_st = 16.7;
    tstart= 2*k-2
    tend=tstart+10*24;
    %y0=[10 10 13000 950 50 0.0333 0 0 0.0333 1];
    %y0=[0.1178 0.1 1.988 0.1553 0.009805 7.5115e-6 0 0 1e-5 1];
    y0=[0.1178 0.1 5 0.6 0.5 7.5115e-6 0 0 1e-5 1];
    tspan=[0,tend];
    
   
    options = odeset('MaxStep',1.0e-2);
    [t, y] = ode15s(@fun_ins, tspan, y0, options);
    % I_0=1.6;    % amol
    I_0=9.2586397e-9;   % ug
    N_c=1000;    N=2.76e6;  % # of pancreatic beta-cells
    F=y(:,6);   G=zeros(length(t),1);       f=zeros(length(t),1);
    toc
    
    for i=1:length(t)
        
        if t(i) < tstart
            G(i)=1;
        elseif mod(t(i)-tstart,24) <= 0
            G(i)=1;
        elseif mod(t(i)-tstart,24) < 1
            G(i)=G_st*0.66667;
        elseif mod(t(i)-tstart,24) <= 5
            G(i)=1;
        elseif mod(t(i)-tstart,24) < 6
            G(i)=G_st*0.66667;
        elseif mod(t(i)-tstart,24) <= 10
            G(i)=1;
        elseif mod(t(i)-tstart,24) < 11
            G(i)=G_st*0.66667;
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
    ISR = I_0 .* sigma .* F .*  f .* N;     % ug/hr
    figure(1)
    vars={'I','V','R','D','D_{IR}','F','\gamma','\rho','F2','G2'};
    for i=1:length(y0)
        hold on
        subplot(4,4,i)
        plot(t,y(:,i))
        xlim([0,tend])
        xlabel('time(h)'); title(vars(i))
        hold on
    end
    
    subplot(4,4,length(y0)+1)
    hold on
    plot(t,ISR)
    xlim([0,tend])
    xlabel('time(h)');    title('ISR');
    
    subplot(4,4,length(y0)+2)
    hold on
    plot(t,G)
    xlim([0,tend])
    xlabel('time(h)');    title('G');
    
    subplot(4,4,length(y0)+3)
    hold on
    plot(t,f)
    xlim([0,tend])
    xlabel('time(h)');    title('f');
    hold on
    
    vars={'I','V','R','D','D_{IR}','F'};
    figure(2)
    for i=1:6
        subplot(2,4,i+1)
        hold on
        plot(t,y(:,i))
        xlim([0,tend])
        xlabel('time(h)'); title(vars(i))
        hold on
    end
    
    subplot(2,4,1)
    hold on
    plot(t,G)
    xlim([0,tend])
    xlabel('time(h)');    title('G');
    
    subplot(2,4,8)
    hold on
    plot(t,ISR)
    xlim([0,tend])
    xlabel('time(h)');    title('ISR');
    
    timetext=num2str(tstart);
    FileName=strcat('equalmeal',timetext);
    save(FileName);
    

    clear

end
toc

FileName='equalmeal';
save(FileName);
FileNameAll=strcat(FileName,'_all');

subplot(2,4,1), hold on
legend('0 hr','2 hr','4 hr','6 hr','8 hr','10 hr','12 hr','14 hr','16 hr','18 hr','20 hr','22 hr')

figure(1)
saveas(gca,FileNameAll);
figure(2)
saveas(gca,FileName);
