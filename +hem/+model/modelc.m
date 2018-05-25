% Control model, light from 6am-6pm, feeding from 6am-6pm 

function dy = modelc(t1, y, s)
global FF tt sigma G_star f_b K_f h_rho G_st tstart


%% Read input
% Intrinsic responses

% Corticosteroid signaling
F = y(1);
mRNA_R_F = y(2);
R_F = y(3); % Glucocorticoid receptor
FR = y(4);
FRn = y(5);

% HPA circadian model 
CRH=y(6);
ACTH=y(7);

% Glucocorticoid receptor in peripheral immune cells (2nd compartment)
F_per=y(8);

% Peripheral clock genes
percry = y(9);
PERCRY = y(10);
nucPERCRY = y(11);
bmal = y(12);
BMAL = y(13);
nucBMAL = y(14);
CLOCKBMAL = y(15);

% Glucocorticoid receptor (Periphery)
R_F_per = y(16);%:229 
FR_per = y(17);%:238
FRn_per = y(18);%:247

% Mineralocorticoid receptor
M_F_per=y(19);
FM_per=y(20);
FMn_per=y(21);
% 
% %Proinflammatory response
mRNA_P = y(22);%:103
P = y(23);%:112
mRNA_R_P = y(24);%:121
R_P = y(25);%:130
PR = y(26);%:139
P_cen = y(27);

%Metabolism
NAD = y(28);
NAM = y(29);
NMN = y(30);
SIRT1 = y(31);
CLOCKBMALSIRT1 = y(32);
NAMPT = y(33);
feed2 = y(34);
feed3 = y(35);
EntF = y(36);
pgc1a = y(37);
PGC1a = y(38);
PGC1aN = y(39);
PGC1aNa = y(40);
FOXO1 = y(41);
Gluc = y(42);
I=y(43); V=y(44);     R=y(45);         D=y(46);     D_IR=y(47);
Fins=y(48); gamma=y(49); rho=y(50);       F2=y(51);    G2=y(52);



% Light for the occasion of HPA entrainment regular step function
if t1<720
    if (mod(t1+s.lst3,24.0000) > 12)
        light=1;
    else
        light=0;
    end
else
    if (mod(t1+s.lst3,24.0000) > 12)
        light=1;
    else
        light=0;
    end
end
% % Feeding 3 meals per day
% if (mod(t1-6,24) > 12)
%     feed=0;
% elseif (mod(t1-8,24) > 12)
%     feed=2;
% elseif (mod(t1-11,24) > 12)
%     feed=0;
% elseif (mod(t1-13,24) > 12)
%     feed=2;
% elseif (mod(t1-16,24) > 12)
%     feed=0;
% elseif(mod(t1-18,24) > 12)
%     feed=2;
% else
%     feed=0;
% end

%Feeding as a regular step function
if t1<528
    if (mod(t1-6,24) < 12)
        feed=1;
    else
        feed=0;
    end
else
    if (mod(t1-6,24) < 12)
        feed=1;
    else
        feed=0;
    end
end
if t1 < tstart
    G=1;
elseif mod(t1-tstart,24) <= 0
    G=1;
elseif mod(t1-tstart,24) < 1
    G=16.7;
elseif mod(t1-tstart,24) <= 5
    G=1;
elseif mod(t1-tstart,24) < 6
    G=16.7;
elseif mod(t1-tstart,24) <= 10
    G=1;
elseif mod(t1-tstart,24) < 11
    G=16.7;
elseif mod(t1-tstart,24) >= 11
    G=1;
end

% if t1 < tstart
%     feed=0; G=1;
% elseif mod(t1-tstart,24) <= 0
%     feed=0; G=1;
% elseif mod(t1-tstart,24) < 1
%     feed=4; G=16.7;
% elseif mod(t1-tstart,24) <= 5
%     feed=0; G=1;
% elseif mod(t1-tstart,24) < 6
%     feed=4; G=16.7;
% elseif mod(t1-tstart,24) <= 10
%     feed=0; G=1;
% elseif mod(t1-tstart,24) < 11
%     feed=4; G=16.7;
% elseif mod(t1-tstart,24) >= 11
%     feed=0; G=1;
% end


% % Light for the occasion of HPA entrainment regular step function
% if (mod(t1+s.lst3,24.0000) > s.sp)
%     light=s.lst1;
% else
%     light=s.lst2;
% end
% % Feeding as a regular step function
% if (mod(t1-s.fdphase,24) < s.fdpd)
%     feed=s.fdhigh;
% else
%     feed=s.fdlow;
% end

%s.new(s.paramtest) = s.new(s.paramtest)*1.01;
dCRH=(s.new(1))/(s.new(2)+FRn)-s.new(3)*CRH*(1+(1*light)/(1+light))/(s.new(4)+CRH) ; % the Hill function of P... *(1+0*s.hh_opt(55)*P^s.hh_opt(57)/(s.hh_opt(56)+P^s.hh_opt(57)))
dACTH=s.new(5)*(1+s.new(22)*P_cen)*CRH/(s.new(2)+FRn)-s.new(6)*ACTH/(s.new(7)+ACTH);
%dACTH=0.4*(1+0.1*P_cen)*CRH/(s.new(2)+FRn)-1*ACTH/(s.new(7)+ACTH);
%dACTH=0.4*(1+s.new(22)*P_cen)*CRH/(s.new(2)+FRn)-s.new(6)*ACTH/(s.new(7)+ACTH);

% dF=s.new(8)*ACTH*(1+(s.new(22)*P_cen)/(s.new(23)+P_cen))-s.new(9)*F/(s.new(10)+F);
% dF=s.new(8)*((ACTH)/(s.new(27)
% +ACTH))*(1+s.new(22)*P_cen)-s.new(9)*F/(s.new(10)+F); Kamau's code, MM
%dF=s.new(8)*ACTH*(1+s.new(22)*P_cen)-s.new(9)*F/(s.new(10)+F); % Pantelis + Pro-inflammatory cytokine
%dF=s.new(8)*ACTH*(1+s.new(22)*P_cen)*(1+0.1*NAD)-s.new(9)*F/(s.new(10)+F); %tryign indirect response
dmRNA_R_F = s.new(11)*(1 - (FRn/(s.new(12)+FRn))) - s.new(13)*mRNA_R_F;
dR_F = s.new(14)*mRNA_R_F + s.new(15)*s.new(16)*FRn - s.new(17)*(F)*R_F - s.new(18)*R_F;
dFR = s.new(17)*(F)*R_F - s.new(19)*FR; % There should be a q in the last term, but then FR just accumulates since there is no degradation/dissociation in this model
dFRn = s.new(19)*FR - s.new(16)*FRn;
% dFRn = 0.83032*FR - s.new(20)*FRn;



% Peripheral glucocorticoid receptor pharmacodynamics (2nd compartment_Clock)
dF_per=(1/(s.new(24)))*(F-F_per);
dP_cen = (1/(s.new(25)))*(P-P_cen);


% Peripheral clock genes
f =1;    % Introduced in order to play with params *random('Uniform',0.8,1.2);

%s.hh_all(s.paramtest) = s.hh_all(s.paramtest)*1.01;

v1b = s.hh_all(41)*f;                % originally is 9 but replaced by 4 in order to simulate human data
k1b = s.hh_all(42)*f;
k1i = s.hh_all(43)*f;              %0.56;
k1d = s.hh_all(44)*f;              %0.18;
k2b = s.hh_all(45)*f;
k2d = s.hh_all(46)*f;              %0.1;
%k3d = s.hh_all(47)*f;              %0.18;
%k2t = s.hh_all(48)*f;              %0.36;
%k3t = s.hh_all(49)*f;
v4b = s.hh_all(50)*f;               %1;
k4b = s.hh_all(51)*f*0.8;
k4d = s.hh_all(52)*f;              %1.1;
k5b = s.hh_all(53)*f;
k5d = s.hh_all(54)*f;              %0.09;
k6d = s.hh_all(55)*f;              %0.18;
k5t = s.hh_all(56)*f;
k6t = s.hh_all(57)*f;              
k6a = s.hh_all(58)*f;               %8
k7a = s.hh_all(59)*f;             %16  
k7d = s.hh_all(60)*f;             %0.13;
c = s.hh_all(61)*f;
p = s.hh_all(62)*f;                   %3 
q = s.hh_all(63)*f;
r = s.hh_all(64)*f;
ent=s.hh_all(65)*f;
% Peripheral clock genes
    
% dpercry = v1b*(1+s.new(26)*P)*(CLOCKBMAL + c)/(k1b*(1+(nucPERCRY/k1i)^p))-...
%     k1d*percry+ ent*FRn_per/(CLOCKBMAL);            
% dPERCRY = k2b*percry^q - k2d*PERCRY - k2t*PERCRY + k3t*nucPERCRY;
%dnucPERCRY = k2t*PERCRY - k3t*nucPERCRY - k3d*nucPERCRY;
%dbmal = (v4b*nucPERCRY^r)/((k4b^r + nucPERCRY^r)) - k4d*bmal;
dBMAL = k5b*bmal - k5d*BMAL - k5t*BMAL + k6t*nucBMAL;
dnucBMAL = k5t*BMAL - k6t*nucBMAL - k6d*nucBMAL + k7a*CLOCKBMAL - k6a*nucBMAL;
%dCLOCKBMAL = k6a*nucBMAL - k7a*CLOCKBMAL - k7d*CLOCKBMAL;


%Receptor Dynamics

%dGR
dR_F_per = s.hh_all(69)*(1+s.hh_all(82)*F_per/(s.hh_all(83)+F_per))*...
    (s.hh_all(70)-R_F_per)/(s.hh_all(71)+s.hh_all(70)-R_F_per)-...
    s.hh_all(72)*R_F_per/(s.hh_all(73)+R_F_per)-...
    s.new(28)*R_F_per*F_per+s.new(29)*FRn_per;
%dFGR
dFR_per = s.hh_all(74)*(F_per)*R_F_per - s.hh_all(75)*FR_per; 
%dFGR(N)
dFRn_per = s.hh_all(75)*FR_per - s.hh_all(76)*FRn_per;

%dMR
dM_F_per=s.hh_all(77)*(1+s.hh_all(84)*F_per/(s.hh_all(85)+F_per))*...
    (s.hh_all(78)-M_F_per)/(s.hh_all(79)+s.hh_all(78)-M_F_per)-...
    s.hh_all(80)*M_F_per/(s.hh_all(81)+M_F_per)-...
    s.new(28)*M_F_per*F_per+s.new(29)*FMn_per;
%dFMR
dFM_per = s.hh_all(74)*(F_per)*M_F_per - s.hh_all(75)*FM_per;
%dFMR(N)
dFMn_per = s.hh_all(75)*FM_per - s.hh_all(76)*FMn_per;

% %Inflammatory Reponse
dmRNA_P = s.hh_all(16)*(1 - (s.hh_all(18)*(FRn_per)/(s.new(30)+ FRn_per)))*(1+PR)*...
    (1-s.hh_all(19)*nucBMAL/(s.new(23)+nucBMAL)) - mRNA_P*s.hh_all(20); 
dP =  mRNA_P*s.hh_all(21) - s.hh_all(22)*P;
dmRNA_R_P = s.hh_all(23)*(1 + (s.hh_all(24)*(FMn_per)/((s.new(31) + FMn_per)))) - s.hh_all(25)*mRNA_R_P;  
dR_P = s.hh_all(26)*mRNA_R_P-s.hh_all(27)*P*R_P-s.hh_all(28)*R_P;  
dPR = s.hh_all(27)*P*R_P - s.hh_all(29)*PR;

% Additional equations for metabolism
%s.met(s.paramtest) = s.met(s.paramtest)*1.01;
nad = s.met(1);	km1 = s.met(2);	Km1 = s.met(3);	km2 = s.met(4);	Km2 = s.met(5);   
km3 = s.met(6);	Km3 = s.met(7);	km4 = s.met(8); Km4 = s.met(9); km5 = s.met(10);
Km5 = s.met(11);    sirtT = s.met(12);  km6 = s.met(13);    Km6 = s.met(14);    km7 = s.met(15);
Km7 = s.met(16);    k3d = s.met(17);    km8a = s.met(18);   km8d = s.met(19);   km9d = s.met(20);
km10a = s.met(21);  km10d = s.met(22);  tau = s.met(23);    k2t = s.met(24);    k3t = s.met(25);    
k6a = s.met(26);    k7a = s.met(27);    k7d = s.met(28);    kn = s.met(29);     km11 = s.met(30); 
Km11 = s.met(31);   km12 = s.met(32);

% KO experiment


dfeed2=1/tau*(feed-feed2);  dfeed3=1/tau*(feed2-feed3);
dNAD = km1*(nad-NAD)/(Km1+nad-NAD) + km2*NMN/(Km2+NMN) - km3*feed3*NAD/(Km3+NAD) - km4*NAD/(Km4+NAD);
dNAM = km4*NAD/(Km4+NAD) - km5*NAMPT*NAM/(Km5+NAM);
dNMN = km5*NAMPT*NAM/(Km5+NAM) - km2*NMN/(Km2+NMN);
dSIRT1 = km6*NAD*(sirtT-SIRT1)/(Km6+sirtT-SIRT1) - km7*SIRT1/(Km7+SIRT1)...
    -km8a*CLOCKBMAL*SIRT1+km8d*CLOCKBMALSIRT1;
dnucPERCRY = k2t*PERCRY - k3t*nucPERCRY - k3d*nucPERCRY*(1+SIRT1);
dCLOCKBMAL = k6a*nucBMAL - k7a*CLOCKBMAL - k7d*CLOCKBMAL...
    -km8a*CLOCKBMAL*SIRT1+km8d*CLOCKBMALSIRT1;
dCLOCKBMALSIRT1 = km8a*CLOCKBMAL*SIRT1 - km8d*CLOCKBMALSIRT1 - km9d*CLOCKBMALSIRT1;
dNAMPT=km10a*CLOCKBMALSIRT1 - km10d*NAMPT;
dpercry = v1b*(1+s.new(26)*P)*(CLOCKBMAL + c)/(k1b*(1+(nucPERCRY/k1i)^p))-...
     k1d*percry+ ent*FRn_per/(CLOCKBMAL);     
%dpercry = 4*(1+0.1*P)*(CLOCKBMAL + c)/(k1b*(1+(nucPERCRY/k1i)^p))-...
%    k1d*percry+ ent*FRn_per/(CLOCKBMAL); 
dPERCRY = k2b*percry^q - k2d*PERCRY - k2t*PERCRY + k3t*nucPERCRY;
dF=s.new(8)*ACTH*(1+s.new(22)*P_cen)+kn*(1+1*EntF/(1+EntF))-s.new(9)*F/(s.new(10)+F);
%dF=0.15*ACTH*(1+s.new(22)*P_cen)+1.5*(1+1*EntF/(1+EntF))-3.4*F/(s.new(10)+F);
dEntF = km11*(NAD)/(Km11+NAD) - km12*EntF;


%s.glu(s.paramtest) = s.glu(s.paramtest)*1.01;
kg1=s.glu(1);  kg2=s.glu(2);  kg3b=s.glu(3); kg3d=s.glu(4); kg3t=s.glu(5);
kg4t=s.glu(6); kg5=s.glu(7);  kg7=s.glu(8);  kg8=s.glu(9);  kg8d=s.glu(10);
kg9=s.glu(11);  kg10=s.glu(12); kg11=s.glu(13); kg12=s.glu(14); kg13=s.glu(15);   
s=s.glu(16);


% kg1=2;  kg2=2;  kg3b=3; kg3d=3; kg3t=2;
% kg4t=2; kg5=0.1;  kg6=1;  kg7=1;  kg8=0.8;
% kg8d=0.5;   kg9=3;  kg10=5; kg11=70; kg12=3; 
% kg13=0.1;   s=8;

dpgc1a = kg1*(1+FOXO1) - kg2*pgc1a;
dPGC1a = kg3b*pgc1a - kg3d*PGC1a - kg3t*PGC1a + kg4t*PGC1aN;
dPGC1aN = kg3t*PGC1a - kg4t*PGC1aN + kg8d*PGC1aNa - kg5*PGC1aN*(SIRT1); %%
dPGC1aNa = kg5*PGC1aN*(1+SIRT1) - kg8*PGC1aNa - kg8d*PGC1aNa;
dFOXO1 = kg9*(PGC1aNa) - kg10*FOXO1;
%dGluc = kg11*FOXO1*(1+kg7*FRn_per) - kg12*Gluc;
dGluc = kg11*FOXO1*FRn_per/(1+(nucPERCRY/kg7)^s) - kg12*Gluc;
dbmal = (0.5*v4b*nucPERCRY^r)/((k4b^r + nucPERCRY^r))*(1+kg13*PGC1aNa) - k4d*bmal;



k=1e-2*60;          alpha_I=0.3*60;     %b_I=4*60;           
tau_G=60/60;
alpha_V=0.6*60;
b_V=bmal;  % med level
sigma=30*60;        tau_V=5/60;
k_1p=5.788e-5*60;   C_T=300;            k_1n=0.255*60;
eta=4*60;           G_star=4.58;        G_hat=10;           h_hat=3.93e-3*60;
rho_b=0.02*60;      gamma_b=1e-4*60;    k_rho=350;          zeta=4*60;
f_b=0.05;           K_f = 3.43;     

dG2 = 1/tau_G*(G-G2);

b_I=60*(2.3.*G2.^3./(4.^3+G2.^3));
dI = -k*I*V - alpha_I*I + b_I;


dV = -k*I*V - alpha_V*V + 10*b_V + sigma*F2;


dR = k*I*V - gamma*R;

dD = gamma*R - k_1p*(C_T - D_IR)*D + k_1n*D_IR;

dD_IR = k_1p*(C_T - D_IR)*D - k_1n*D_IR - rho*D_IR;

dFins = rho*D_IR - sigma*Fins;
dF2 = 1/tau_V*(Fins-F2);

if G2 <= G_star
    h_gamma = 0;
elseif G2 <= G_hat
    h_gamma = h_hat*(G2-G_star)/(G_hat-G_star);
else
    h_gamma = h_hat;
end

dgamma = eta*(-gamma + gamma_b + h_gamma);

if gamma < gamma_b
    h_rho = 0;
else
    h_rho = k_rho*(gamma - gamma_b);
end

drho = zeta*(-rho + rho_b + h_rho);

%% Format output

dy = zeros(size(y));

% Intrinsic responses

% Corticosteroid signaling
dy(1) = dF;
dy(2) = dmRNA_R_F;
dy(3) = dR_F;
dy(4) = dFR;
dy(5) = dFRn;


% HPA circadian model
dy(6)=dCRH;
dy(7)=dACTH;

% Peripheral Cortisol
dy(8)=dF_per;

% % Peripheral clock genes
dy(9) = dpercry;
dy(10) = dPERCRY;
dy(11) = dnucPERCRY;
dy(12) = dbmal;
dy(13) = dBMAL;
dy(14) = dnucBMAL;
dy(15) = dCLOCKBMAL;

%Receptor Dynamics
dy(16) = dR_F_per;
dy(17) = dFR_per;
dy(18) = dFRn_per;
dy(19) = dM_F_per;
dy(20) = dFM_per;
dy(21) = dFMn_per;
% 
% %Inflammatory Module
dy(22) = dmRNA_P;
dy(23) = dP;
dy(24) = dmRNA_R_P;
dy(25) = dR_P;
dy(26) = dPR;
dy(27) = dP_cen;

% Metabolism
dy(28) = dNAD;
dy(29) = dNAM;
dy(30) = dNMN;
dy(31) = dSIRT1;
dy(32) = dCLOCKBMALSIRT1;
dy(33) = dNAMPT;
dy(34) = dfeed2;
dy(35) = dfeed3;
dy(36) = dEntF;
dy(37) = dpgc1a;
dy(38) = dPGC1a;
dy(39) = dPGC1aN;
dy(40) = dPGC1aNa;
dy(41) = dFOXO1;
dy(42) = dGluc;

dy(43)=dI;   dy(44)=dV;       dy(45)=dR;       dy(46)=dD;   dy(47)=dD_IR;
dy(48)=dFins;   dy(49)=dgamma;   dy(50)=drho;     dy(51)=dF2;  dy(52)=dG2;
end