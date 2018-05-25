%cortisol model
 function dy = cortisol_only(t1, y, s)
    F = y(1);
    mRNA_R_F = y(2);
    R_F = y(3); % Glucocorticoid receptor
    FR = y(4);
    FRn = y(5);
    CRH = y(6);
    ACTH = y(7);
    
    light = 0;
    EntF = 0;
    P_cen = 0;
    kn=0.7;
    
dCRH=(s.new(1))/(s.new(2)+FRn)-s.new(3)*CRH*(1+(1*light)/(1+light))/(s.new(4)+CRH) ; % the Hill function of P... *(1+0*s.hh_opt(55)*P^s.hh_opt(57)/(s.hh_opt(56)+P^s.hh_opt(57)))
dACTH=s.new(5)*(1+s.new(22)*P_cen)*CRH/(s.new(2)+FRn)-s.new(6)*ACTH/(s.new(7)+ACTH);
dmRNA_R_F = s.new(11)*(1 - (FRn/(s.new(12)+FRn))) - s.new(13)*mRNA_R_F;
dR_F = s.new(14)*mRNA_R_F + s.new(15)*s.new(16)*FRn - s.new(17)*(F)*R_F - s.new(18)*R_F;
dFR = s.new(17)*(F)*R_F - s.new(19)*FR; % There should be a q in the last term, but then FR just accumulates since there is no degradation/dissociation in this model
dFRn = s.new(19)*FR - s.new(20)*FRn;
%dF=s.new(8)*ACTH*(1+s.new(22)*P_cen)*(1+kn*EntF)-s.new(9)*F/(s.new(10)+F);
dF=0.5*ACTH*(1+s.new(22)*P_cen)*(1+kn*EntF)-s.new(9)*F/(s.new(10)+F);
dy = zeros(size(y));
dy(1) = dF;
dy(2) = dmRNA_R_F;
dy(3) = dR_F;
dy(4) = dFR;
dy(5) = dFRn;


% HPA circadian model
dy(6)=dCRH;
dy(7)=dACTH;

 end