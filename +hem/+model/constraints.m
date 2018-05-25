% clear all;% close all;
% clc;
% tic;

function [c,ceq]=constraints(hh_change,hh_store,v, s)

for cu=1:length(v)
    hh_store(v(cu))=hh_change(cu);
end
hh_store(64)
hh_store(65)
s.hh_opt = hh_store;
s.make_plots=false;

[t, y] = hem.model.run(s);
    
for tim=1:length(t)
        if (mod(t(tim)+s.lst3,24) >s.sp)
            light(tim)=hh_store(61);
        else
            light(tim)=s.lst2;
        end
    end
    %% Calculate period/phase with FFT
    % Prepare the data for FFT analysis.
    Fs = 10;    % Sampling frequency that will be used in the FFT.
    b = (0:1/Fs:(s.t_final))+s.t_pre;   % Interpolation of signal in order to be used in FFT.
    e1 = interp1(t,y.F_clock,b);
    e2=interp1(t,y.percry,b);
    e3=interp1(t,light,b);
    e10 = interp1(t,y.P,b);
    e11 = interp1(t,y.A,b);
    
    mv1=mean(e1);
    mv2=mean(e2);
    mv3=mean(e3);

    mx=max(e1);
    mn=min(e1);
    mx2=max(e2);
    mn2=min(e2);

    E1 = fft(e1);
    [zud1,m1]=max(abs(E1(2:end)));
    ang11 = angle((E1));
    ang1=ang11(m1);
    while ang1<0
        ang1=ang1+2*pi;
    end
    
     E2 = fft(e2);
    [zud2,m2]=max(abs(E2(2:end)));
    ang22 = angle((E2));
    ang2=ang22(m2);
    while ang2<0
        ang2=ang2+2*pi;
    end
        
    E3 = fft(e3);
    [zud3,m3]=max(abs(E3(2:end)));
    ang33 = angle((E3));
    ang3=ang33(m3);
    while ang3<0
        ang3=ang3+2*pi;
    end
    
    T1=(m1*Fs/2)/(s.t_final*Fs/2)
    T2=(m2*Fs/2)/(s.t_final*Fs/2)
    T3=(m3*Fs/2)/(s.t_final*Fs/2)

    
    c(1) = -hh_store(64);
    c(2) = -hh_store(65);
    ceq = T1-0.0408-hh_store(64)+hh_store(65);
end

