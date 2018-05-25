function dy = fun_ins(t1,y)

global FF tt sigma G_star f_b K_f h_rho G_st tstart

I=y(1); V=y(2);     R=y(3);         D=y(4);     D_IR=y(5);
F=y(6); gamma=y(7); rho=y(8);       F2=y(9);    G2=y(10);

k=1e-2*60;          alpha_I=0.3*60;     %b_I=4*60;           
tau_G=60/60;
alpha_V=0.6*60;     %b_V=6*60;           
%b_V=3*sin((2*3.14/24)*t1-(2*3.14/24*12))+6;  %original
%b_V=5*sin((2*3.14/24)*t1-(2*3.14/24*12))+6;  %high amplitude
%b_V=3*sin((2*3.14/24)*t1-(2*3.14/24*12))+100;  %level up
%b_V=5*sin((3.14/12)*(t1-10))+6;  %high amplitude (all results prior to
%1/18/18
b_V=3*sin((3.14/12)*(t1-10))+9;  % med level
sigma=30*60;        tau_V=5/60;
k_1p=5.788e-5*60;   C_T=300;            k_1n=0.255*60;
eta=4*60;           G_star=4.58;        G_hat=10;           h_hat=3.93e-3*60;
rho_b=0.02*60;      gamma_b=1e-4*60;    k_rho=350;          zeta=4*60;
f_b=0.05;           K_f = 3.43;     




if t1 < tstart
    G=1;
elseif mod(t1-tstart,24) <= 0
    G=1;
elseif mod(t1-tstart,24) < 1
    G=G_st*0.66667;
elseif mod(t1-tstart,24) <= 5
    G=1;
elseif mod(t1-tstart,24) < 6
    G=G_st*0.66667;
elseif mod(t1-tstart,24) <= 10
    G=1;
elseif mod(t1-tstart,24) < 11
    G=G_st*0.66667;
elseif mod(t1-tstart,24) >= 11
    G=1;
end

% if (mod(t1-6,24.0000) < 12)
%     G=G_st;
% else
%     G=1;
% end
% 
% if (t1 <= 2+24)
%     G=1;
% elseif (t1 < 3+24)
%     G=G_st;
% elseif (t1 <= 5+24)
%     G=1;
% elseif (t1 < 6+24)
%     G=G_st;
% elseif (t1 <= 8+24)
%     G=1;
% elseif (t1 < 9+24)
%     G=G_st;
% elseif (t1 <= 11+24)
%     G=1;
% elseif (t1 < 12+24)
%     G=G_st;
% elseif (t1 <= 14+24)
%     G=1;
% elseif (t1 < 15+24)
%     G=G_st;
% elseif (t1 <= 17+24)
%     G=1;
% elseif (t1 < 18+24)
%     G=G_st;
% elseif (t1 <= 20+24)
%     G=1;
% elseif (t1 < 21+24)
%     G=G_st;
% elseif (t1 <= 23+24)
%     G=1;
% elseif (t1 < 24+24)
%     G=G_st;
% else
%     G=1;
% end

dG2 = 1/tau_G*(G-G2);

b_I=60*(2.3.*G2.^3./(4.^3+G2.^3));
dI = -k*I*V - alpha_I*I + b_I;


dV = -k*I*V - alpha_V*V + b_V + sigma*F2;


dR = k*I*V - gamma*R;

dD = gamma*R - k_1p*(C_T - D_IR)*D + k_1n*D_IR;

dD_IR = k_1p*(C_T - D_IR)*D - k_1n*D_IR - rho*D_IR;

dF = rho*D_IR - sigma*F;
dF2 = 1/tau_V*(F-F2);

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


dy = zeros(size(y));
dy(1)=dI;   dy(2)=dV;       dy(3)=dR;       dy(4)=dD;   dy(5)=dD_IR;
dy(6)=dF;   dy(7)=dgamma;   dy(8)=drho;     dy(9)=dF2;  dy(10)=dG2;

FF=[FF F];
tt=[tt t1];
end