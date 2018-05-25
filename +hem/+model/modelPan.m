 function dy = model(t, y, s, kt, n)
%HEM.MODEL.MODEL Evaluate human endotoxemia model equations.
% Note: Although Proinflammatory(P) and anti-inflammatory responsevariables
% are defined they are not used in the HPA AXIS equations
%
%_________________________________________
% s.hh(10,1) =s.hh(8,1)*(1+s.hh(9,1));
% s.hh(15,1) = s.hh(13,1);
% s.hh(21,1) = s.hh(18,1)*(1+s.hh(20,1));
% s.hh_lps(4,1) = s.hh_lps(1,1)*(1+s.hh_lps(2,1))*(1+s.hh_lps(3,1));
% s.hh(28,1) = s.hh(26,1)*(1+s.hh(27,1));
% s.hh_all(3,1) = s.hh_all(1,1)*(1+s.hh_all(2,1));
% s.hh_all(5,1) = s.hh_all(6,1)*(1+s.hh_all(7,1)) + s.hh_all(8,1);
% s.hh_fen(3,1) = s.hh_fen(1,1)*(1+s.hh_fen(2,1));
%_________________________________________________


%% Read input
% Intrinsic responses
P = y(1); % Pro-inflammatory response
A = y(2); % Anti-inflammatory response

% Corticosteroid signaling
F = y(3);
mRNA_R_F = y(4);
R_F = y(5); % Glucocorticoid receptor
FR = y(6);
FRn = y(7);

% HPA circadian model 
CRH=y(8);
ACTH=y(9);

% Glucocorticoid receptor in peripheral immune cells (2nd compartment)
F_clock=y(10);
mRNA_R_F_clock=y(11);
R_F_clock=y(12);
FR_clock=y(13);
FRn_clock=y(14);

% Peripheral clock genes
percry = y(15);
PERCRY = y(16);
nucPERCRY = y(17);
bmal = y(18);
BMAL = y(19);
nucBMAL = y(20);
CLOCKBMAL = y(21);

% Pro-Inflammatory cytokines in the central clock (1st compartment)
P_central=y(22);

%% Differential equations

% Intrinsic responses

dP = s.hh_opt(14)*(1+s.hh_opt(54)*percry)/(A) - s.hh_opt(17)*P;

dA = s.hh_opt(18)*(1+s.hh_opt(20)*FRn_clock/(CLOCKBMAL))*(1+s.hh_opt(21)*percry) - s.hh_opt(22)*A;

% Light for the occasion of HPA entrainment regular step function
% if (mod(t+s.lst3,24.0000) > s.sp)
%     light=s.lst1;
% else
%     light=s.lst2;
% end
% 
% if (mod(t+kt,24.0000) > 22)
%        light2=1*n;
% else
%        light2=0;
% end

% if (mod(t,24.000) < 2)
%     light=0;
% elseif (mod(t,24.0000) < 3)
%     light=1;
% else
if (mod(t,24.0000) < 6)
    light=0;
elseif (mod(t,24.000) < 18)
    light=1;
elseif (mod(t,24.000) < 19)
    light=0;
elseif (mod(t,24.000) < 21)
    light=1;
else
    light=0;
end

% Less parameters without linear terms
dCRH=s.hh_opt(1)/(s.hh_opt(3)+FRn)-s.hh_opt(4)*CRH*(1+s.hh_opt(61)*(light)/(1+(light)))/(s.hh_opt(10)+CRH);  % the Hill function of P... *(1+0*s.hh_opt(55)*P^s.hh_opt(57)/(s.hh_opt(56)+P^s.hh_opt(57)))
dACTH=s.hh_opt(2)*CRH/(s.hh_opt(3)+FRn)-s.hh_opt(6)*ACTH/(s.hh_opt(11)+ACTH);
dF=s.hh_opt(7)*ACTH-s.hh_opt(8)*F/(s.hh_opt(12)+F);

% Glucocorticoid receptor pharmacodynamics (1st compartment_HPA)
dmRNA_R_F = s.hh_jusko(1,1)*(1 - (FRn/(s.hh_jusko(2,1)+FRn))) - s.hh_jusko(3,1)*mRNA_R_F;
dR_F = s.hh_jusko(4,1)*mRNA_R_F + s.hh_jusko(5,1)*s.hh_jusko(6,1)*FRn - s.hh_jusko(7,1)*(F)*R_F - s.hh_jusko(12,1)*R_F;
dFR = s.hh_jusko(7,1)*(F)*R_F - s.hh_jusko(13,1)*FR; 
dFRn = s.hh_jusko(13,1)*FR - s.hh_jusko(14,1)*FRn;


% Peripheral glucocorticoid receptor pharmacodynamics (2nd compartment_Clock)
dF_clock=(1/(s.hh_opt(26)))*(F-F_clock);
dmRNA_R_F_clock = s.hh_jusko(1,1)*(1 - (FRn_clock/(s.hh_jusko(2,1)+FRn_clock))) - s.hh_jusko(3,1)*mRNA_R_F_clock; % 2 CLOCKBMAL1 missing
dR_F_clock = s.hh_jusko(4,1)*mRNA_R_F_clock + s.hh_jusko(5,1)*s.hh_jusko(6,1)*FRn_clock - s.hh_jusko(7,1)*(F_clock)*R_F_clock - s.hh_jusko(12,1)*R_F_clock;
dFR_clock = s.hh_jusko(7,1)*(F_clock)*R_F - s.hh_jusko(13,1)*FR_clock; % There should be a q in the last term, but then FR just accumulates since there is no degradation/dissociation in this model
dFRn_clock = s.hh_jusko(13,1)*FR_clock - s.hh_jusko(14,1)*FRn_clock;



% Peripheral clock genes
f = 1;    % Introduced in order to play with params *random('Uniform',0.8,1.2);
v1b = s.hh_opt(29)*f;                % originally is 9 but replaced by 4 in order to simulate human data
k1b = s.hh_opt(30)*f;
k1i = s.hh_opt(31)*f;              %0.56;
k1d = s.hh_opt(32)*f;              %0.18;
k2b = s.hh_opt(33)*f;
k2d = s.hh_opt(34)*f;              %0.1;
k3d = s.hh_opt(35)*f;              %0.18;
k2t = s.hh_opt(36)*f;              %0.36;
k3t = s.hh_opt(37)*f;
v4b = s.hh_opt(38)*f;               %1;
k4b = s.hh_opt(39)*f;
k4d = s.hh_opt(40)*f;              %1.1;
k5b = s.hh_opt(41)*f;
k5d = s.hh_opt(42)*f;              %0.09;
k6d = s.hh_opt(43)*f;              %0.18;
k5t = s.hh_opt(44)*f;
k6t = s.hh_opt(45)*f;              
k6p = s.hh_opt(46)*f;               %8
k7p = s.hh_opt(47)*f;             %16  
k7d = s.hh_opt(48)*f;             %0.13;
c = s.hh_opt(49)*f;
p = s.hh_opt(50)*f;                   %3 
q = s.hh_opt(51)*f;
r = s.hh_opt(52)*f;
ent=s.hh_opt(53)*f;
  
dpercry = v1b*(CLOCKBMAL + c)/(k1b*(1+(nucPERCRY/k1i)^p))- k1d*percry+ ent*FRn_clock/(CLOCKBMAL) ;            
dPERCRY = k2b*percry^q - k2d*PERCRY - k2t*PERCRY + k3t*nucPERCRY;
dnucPERCRY = k2t*PERCRY - k3t*nucPERCRY - k3d*nucPERCRY;
dbmal = (v4b*nucPERCRY^r)/((k4b^r + nucPERCRY^r)) - k4d*bmal ;
dBMAL = k5b*bmal - k5d*BMAL - k5t*BMAL + k6t*nucBMAL;
dnucBMAL = k5t*BMAL - k6t*nucBMAL - k6d*nucBMAL + k7p*CLOCKBMAL - k6p*nucBMAL;
dCLOCKBMAL = k6p*nucBMAL - k7p*CLOCKBMAL - k7d*CLOCKBMAL;

% Pro-Inflammatory cytokines in the central compartment
dP_central=(1/(s.hh_opt(26)))*(P-P_central);

%% Format output

dy = zeros(size(y));

% Intrinsic responses
dy(1) = dP;
dy(2) = dA;

% Corticosteroid signaling
dy(3) = dF;
dy(4) = dmRNA_R_F;
dy(5) = dR_F;
dy(6) = dFR;
dy(7) = dFRn;

% HPA circadian model
dy(8)=dCRH;
dy(9)=dACTH;

% Peripheral glucocorticoid receptor pharmacodynamics
dy(10)=dF_clock;
dy(11)=dmRNA_R_F_clock;
dy(12)=dR_F_clock;
dy(13)=dFR_clock;
dy(14)=dFRn_clock;

% Peripheral clock genes
dy(15) = dpercry;
dy(16) = dPERCRY;
dy(17) = dnucPERCRY;
dy(18) = dbmal;
dy(19) = dBMAL;
dy(20) = dnucBMAL;
dy(21) = dCLOCKBMAL;

% Pro-Inflammatory cytokines in the central compartment
dy(22)=dP_central;

