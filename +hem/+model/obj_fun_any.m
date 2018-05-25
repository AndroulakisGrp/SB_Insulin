function err = obj_fun_any(hh_change,hh_store,vars,v,s)
% hh_circ = hh_circ_0;
% hh_circ(v) = k



for cu=1:length(v)
    hh_store(v(cu))=hh_change(cu);
end

hh_store
s.hh_all = hh_store;
% keyboard;

figure(1);
clf;


[t, y] = hem.model.run2(s);
Fs = 10;
b = (0:1/Fs:(s.t_final))+s.t_pre; 
%mrnapi = interp1(t,y.mRNA_P,b);
fgi = interp1(t,y.R_F_per,b);
fmi=interp1(t,y.M_F_per,b);

% % Find the amplitudes and mean values of interpolated profiles
% Amps
ampfgi=max(fgi)-min(fgi);
ampfmi=max(fmi)-min(fmi);
%amppri=max(pri)-min(pri);

% Means
meanfgi=mean(fgi);
meanfmi=mean(fmi);
% meanpri=mean(pri);


error = [];

%% Calculate error
% % Error from experimental data fitting
error01=[];

for kai=1:length(vars)
    var = vars{kai};
    [data_t, data_y] = hem.model.get_data(var);
%     if var== 'F'
%         weight=2;
%     else
    weight=1;
%     end
    for toy=1:length(data_t)
        %     if toy==3
        %         f=50;
        %     else
        %             f=1;
        %         else if toy==6
        %                 f=2;
        %             else
        %                 f=1;
        %             end
        %         end
        %     end
        %     err = [err, f*((interp1(t, [y.(var)], data_t(toy)+s.t_pre)-data_y(toy))^2)];
        % koita(toy,kai)=interp1(t, [y.(var)], data_t(toy)+s.t_pre)
        % des(toy,kai)=data_y(toy)
        error01(kai,toy) = ((interp1(t, [y.(var)], data_t(toy)+s.t_pre)-data_y(toy)))^2;
    end
    error1(kai)=weight*sum(error01(kai,:));
    error01=[];
end




% Error produced by acute stress response
figure(2);
clf;
s.LPS=1;
s.t_LPS=21;
[t, y] = hem.model.run2(s);
for kai=1:length(vars)
    var = vars{kai};
    [data_t, data_y] = hem.model.get_data_2(var);
    weight=1;
    for toy=1:length(data_t)
        %     if toy==3
        %         f=50;
        %     else
        %             f=1;
        %         else if toy==6
        %                 f=2;
        %             else
        %                 f=1;
        %             end
        %         end
        %     end
        %     err = [err, f*((interp1(t, [y.(var)], data_t(toy)+s.t_pre)-data_y(toy))^2)];
        % koita(toy,kai)=interp1(t, [y.(var)], data_t(toy)+s.t_pre)
        % des(toy,kai)=data_y(toy)
        error02(kai,toy) = ((interp1(t, [y.(var)], data_t(toy)+s.t_pre)-data_y(toy)))^2;
    end
    error2(kai)=weight*sum(error02(kai,:));
    error02=[];
end


% [t, y] = hem.model.run(0.1, 9, s);
% err = [];
% for sai=1:length(vars)
%     var = vars{sai};
%     [data_t, data_y] = hem.model.get_data_2(var);
% % end
% if var =='P'
%     weight=3;
% else
%     weight=1;
% end
% 
%     
% %     if ~isempty(data_t)
% %             if strcmp('F', var)
% %                 f = 2;
% %             else
% %                 f = 1;
% %             end
% for soy=1:length(data_t)
% %     if toy==3
% %         f=50;
% %     else  if toy==3
% %             f=2;
% %         else if toy==6
% %                 f=2;
% %             else
% %                 f=1;
% %             end
% %         end
% %     end  
% %     err = [err, f*((interp1(t, [y.(var)], data_t(toy)+s.t_pre)-data_y(toy))^2)];
% % koita(toy,kai)=interp1(t, [y.(var)], data_t(toy)+s.t_pre)
% % des(toy,kai)=data_y(toy)
%     error02(soy,sai) = ((interp1(t, [y.(var)], data_t(soy)+s.t_pre)-data_y(soy)))^2;
% end
% error2(sai)=weight*sum(error02(:,sai));
% end


% % Error produced when P is not present in HPA network
% f= hem.util.get_hem_folder;
% hh_circ_3=importdata(sprintf('%s/param/hh_opt.txt', f));
% hh_circ_3(55)=0;
% hh_circ_3
% s.hh_opt=hh_circ_3;
% [t, y] = hem.model.run(0, 0, s);
%  Fs = 10;    % Sampling frequency that will be used in the FFT.
%     b = (0:1/Fs:(s.t_final))+s.t_pre;   % Interpolation of signal in order to be used in FFT.
%     e1 = interp1(t,y.F,b);
% e1 = e1 - mean(e1);
% NFFT = 20000;
%     E1 = fft(e1,NFFT);
%     m=max(E1);
%     ang = angle(m);
% 
% %     Calculate the frequency that FFT correspond with
%     f = Fs*linspace(0,1,NFFT);
%     %
%     % Calculate the maximum of the fourier transformation and find the time
%     % that this happens.
%     [o1,l1] = max(E1(1:NFFT/2));
%     T01 = f(l1)



% % Error produced in chronic stress
% figure(3);
% clf;
% s.transition_logistic=true;
% s.t_ex_LPS=[0 400];
% s.t_LPS=0;
% [t, y] = hem.model.run(0.1, 9, s);
% Fs = 10;    % Sampling frequency that will be used in the FFT.
%     b = (300:1/Fs:(s.t_final))+s.t_pre;   % Interpolation of signal in order to be used in FFT.
%     e1c = interp1(t,y.F,b);
% 
% % Remove the bias (subtract the mean value).
%     e1c = e1c - mean(e1c);
%     mxc=max(e1c);
%     mnc=min(e1c);  


drawnow;


err = (ampfmi-0.1)^2 + (meanfmi-1.3)^2+sum(error1)


%+5*(1-meanpri)^2+5*(1-amppri)^2


% fileID=fopen('storobopt.txt','a+');
% fprintf(fileID,'%.4f \n',err);
% fclose(fileID);
% fileID=fopen('storparamopt.txt','a+');
% fprintf(fileID,'%6sf\n','hh_opt');
% fprintf(fileID,'%6.4f\n',hh_store);
% fclose(fileID);

if isnan(err)
    err = 10^6;
end
if ~isreal(err)
    err = 10^6;
end
%     else
%         err = 10^6;
%     end
% keyboard
