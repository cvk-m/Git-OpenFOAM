function [X,fit,pfit]=fobj(X,Lb,Ub)
% Correcting the design vector if it is not within the defined range.
for i=1:size(X,2)
if X(i)>Ub(i)
X(i)=Ub(i);
end
if X(i)<Lb(i)
X(i)=Lb(i);
end
end
% Calculate inequality constraints (g(i)). Number of inequalityconstraints(l) is 4.
g(1)=cokrigingpredictor(X)-cokrigingpredictor1(X);
g(2)=cokrigingpredictor2(X);

% Calculate the cost function (fit).
fit=0.8*cokrigingpredictor(X)/cokrigingpredictor1(X)+0.2*cokrigingpredictor2(X);
% Defining the penalty parameter (represented here with nou). Notice thatpenalty parameter is considered as a very big number and equal for all fourinequality constraints.
nou=10^9;
penalty=0;
for i=1:size(g,2)
if g(i)>0
penalty=penalty+nou*g(i);
end
end
% Calculate the penalized cost function (pfit) by adding measure of penalty  function (penalty).
pfit=fit+penalty;