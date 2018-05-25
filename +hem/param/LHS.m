r = 0.05*p;
lb = p-r;
ub = p+r;
X = lhsdesign(1000,length(p));
D = bsxfun(@plus,lb',bsxfun(@times,X,(ub-lb)'));
% new = D';


hhall = D';
%%
save hh_all.txt hhall -ascii
%%
save ICs.txt  -ascii
%%
save new.txt  -ascii