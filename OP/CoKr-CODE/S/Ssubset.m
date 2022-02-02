function [Xs,Xr]=Ssubset(X,ns)
% Given a sampling plan, returns a subset of a given size with
% optimized space – filling properties (as per the Morris–Mitchell
% criterion).
%
% Inputs:
% X – full sampling plan
% ns – size of the desired subset
%
% Outputs:
% Xs – subset with optimized space – filling properties
% Xr – remainder X\Xs
n=size(X,1);
% Norm and quality metric exponent – different values can be used if
% required
p=1; q=5;
r=randperm(n);
Xs=X(r(1:ns),:);
Xr=X(r(ns+1:end),:);
for j=1:ns
orig_crit=Smmphi(Xs,q,p);
orig_point=Xs(j,:);
% We look for the best point to substitute the current one with
bestsub=1;
bestsubcrit=Inf;
for i=1:n-ns
% We replace the current, jth point with each of the
% remaining points, one by one
Xs(j,:)=Xr(i,:);
crit=Smmphi(Xs,q,p);
if crit< bestsubcrit
bestsubcrit=crit;
bestsub=i;
end
end
if bestsubcrit<orig_crit
Xs(j,:)=Xr(bestsub,:);
else
Xs(j,:)=orig_point;
end
end