clear all;
close all;
clc;

%close all;

ampdatacombined=[];
phasedatacombined=[];
tic
for index=1
clear ydata adata y amp_phase_interp acrophase e1-cat e1-spr
s.t_pre=480;
s.t_final=1000;
s.make_plots=true;
s.extra_plots=true;
s.ss_plots=true;
s.stochastic=false;
s.color='k';
s.plot_data=false;
s.nonnegative=false;
s.HPA=true;
s.transition_logistic=false;
s.paramtest = index;
% s.t_ex_LPS=[200 1800];

fdpd=12;     %12
fdhigh=1;   %1
fdlow=0;    %0
fdphase=6;
%fdphase=index*2;  %6

s.fdhigh=fdhigh;
s.fdlow=fdlow;
s.fdphase=fdphase;
s.fdpd = fdpd;




% Generate an ensemble of cells
cells=1;

scp=[12];   %light duration
lst1=1;%[1]; %light high
lst2=0;%[0]; %light low
lst3=[6];   %light phase

color={'k','b','m','r','g','y','k','b','m','r','g','y','k','b','m','r','g','y','k','b','m','r','g','y','k','b','m','r','g','y'};
% color={[0.2 1 1],[0.2 0.6 1],[0.2 0.2 1],[0.8 0.2 1],[1 0.2 0.4],[1 0.2 0.2],[1 0.6 0.2]};
for tu=1:cells
    tu;
    s.lst1=lst1;
    s.lst2=lst2;
    s.lst3=lst3;
    s.sp=scp;
    s.color=color{1};

    [t, y] = hem.model.run(s);

end

%%
%Amplitude and Phase calculation as well as fitting to cosinor model

%Partition Output for each season
l1 = 1: length(t1);

%Loop through variables to determine amplitude, phase, period etc, we would
%save the interpolated values for calculation but also use these as guesses
%for our fit
%
%
%%
fields = fieldnames(y);
amp_phase_interp = struct;
rhy_interp = struct;
model_int = cell(1,4);
cos_model = cell(4,numel(fields));
periodP = ones([numel(fields), 1]);
phaseP = ones([numel(fields) 1]);
ii=1;
for ii = 1:numel(fields)
%     if ii<=32
%         ii
%     end
    int = y.(fields{ii});
    spr = int(l1);

    %Normalize time to first time point
    tnew1 = t1 - t1(1);
    b = (0:1/100:(100))+720;     % create a vector starting at 960, ending at 1060
    e1_spr = interp1(tnew1,spr,b);
    L = length(b);
    
    %Concatenate
    e1_cat = e1_spr;
    maxf = max(e1_cat,[],2);
    minf =min(e1_cat,[],2);
    ampf = maxf - minf;
    jj=1;
    for jj = 1:size(e1_cat,1);
        fill = e1_cat(jj,:);
        [locs, pks]=hem.util.peakfinder(fill);
        nn=1; % indicating the index of each peak in the cell
        kk=1;
        for kk=1:length(pks)
            if pks(kk)>mean(fill)
                peak(nn)=pks(kk);
                inofpeak(nn)=locs(kk);
                nn=nn+1;
            end
        end
        
        periodP(ii)=b(inofpeak(3))-b(inofpeak(2));
        phaseP(ii)=2*pi*1/periodP(ii)*(b(inofpeak(2))-b(1));
        if phaseP(ii) >2*pi
            phaseP(ii)=phaseP(ii)-2*pi;
        end
        acrophase(jj,1) = phaseP(ii);
    end

    
    %save results from interpolation
    
    interp_out = cat(2,ampf, acrophase);
    amp_phase_interp.(fields{ii}) = interp_out;
    rhy_interp.(fields{ii}) = e1_cat;

end
    for mm = 1:numel(fields)
        ampdata(mm) = amp_phase_interp.(fields{mm})(1);
        phasedata(mm) = amp_phase_interp.(fields{mm})(2);
    end
    
    ampdatacombined = [ampdatacombined ampdata'];
    phasedatacombined = [phasedatacombined phasedata'];

index
end
toc

