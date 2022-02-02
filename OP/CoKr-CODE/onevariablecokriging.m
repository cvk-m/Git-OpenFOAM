% One variable co-Kriging example (similar to Figure 8.1)
%
% Copyright 2007 A I J Forrester
%
% This program is free software: you can redistribute it and/or modify  it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or any
% later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License and GNU
% Lesser General Public License along with this program. If not, see
% <http://www.gnu.org/licenses/>.
%% MODELINFO CONSTRUCTION FOR OBJ 0 1 2
%% OBJ 0 MODELINFO
global ModelInfo
% Expensive points
ModelInfo.Xe=[0; 0.4; 0.6;1];
% Cheap points
ModelInfo.Xc=[0.1;0.2;0.3;0.5;0.7;0.8;0.9;0;0.4;0.6;1];
k=1;
% Calulate expensive observations
for i=1:size(ModelInfo.Xe,1)
ModelInfo.ye(i,1)=onevar(ModelInfo.Xe(i));
end
% Calculate cheap observations
for i=1:size(ModelInfo.Xc,1)
ModelInfo.yc(i,1)=cheaponevar(ModelInfo.Xc(i));
end
%% OBJ 1 MODELINFO
global ModelInfo1
% Expensive points
ModelInfo1.Xe=[0; 0.4; 0.6;1];
% Cheap points
ModelInfo1.Xc=[0.1;0.2;0.3;0.5;0.7;0.8;0.9;0;0.4;0.6;1];
k=1;
% Calulate expensive observations
for i=1:size(ModelInfo1.Xe,1)
ModelInfo1.ye(i,1)=onevar(ModelInfo1.Xe(i));
end
% Calculate cheap observations
for i=1:size(ModelInfo1.Xc,1)
ModelInfo1.yc(i,1)=cheaponevar(ModelInfo1.Xc(i));
end
%% OBJ 2 MODELINFO
global ModelInfo2
% Expensive points
ModelInfo2.Xe=[0; 0.4; 0.6;1];
% Cheap points
ModelInfo2.Xc=[0.1;0.2;0.3;0.5;0.7;0.8;0.9;0;0.4;0.6;1];
k=1;
% Calulate expensive observations
for i=1:size(ModelInfo2.Xe,1)
ModelInfo2.ye(i,1)=onevar(ModelInfo2.Xe(i));
end
% Calculate cheap observations
for i=1:size(ModelInfo2.Xc,1)
ModelInfo2.yc(i,1)=cheaponevar(ModelInfo2.Xc(i));
end
%%
%% MODEL PARAMETERS CONSTRUCTION FOR OBJ 0 1 2
%%
%% OBJ 0 MODEL PARAMETER
% Optimize cheap model paramters
ModelInfo.Thetac=fminbnd(@likelihoodc,-3,3);
% Optimise difference model paramters
Params=ga(@likelihoodd,k+1,[],[],[],[],[-3 -5],[3 5]);
ModelInfo.Thetad=Params(1:k);
ModelInfo.rho=Params(k+1);
%% OBJ 1 MODEL PARAMETER
% Optimize cheap model paramters
ModelInfo1.Thetac=fminbnd(@likelihoodc1,-3,3);
% Optimise difference model paramters
Params=ga(@likelihoodd1,k+1,[],[],[],[],[-3 -5],[3 5]);
ModelInfo1.Thetad=Params(1:k);
ModelInfo1.rho=Params(k+1);
%% OBJ 2 MODEL PARAMETER
% Optimize cheap model paramters
ModelInfo2.Thetac=fminbnd(@likelihoodc2,-3,3);
% Optimise difference model paramters
Params=ga(@likelihoodd2,k+1,[],[],[],[],[-3 -5],[3 5]);
ModelInfo2.Thetad=Params(1:k);
ModelInfo2.rho=Params(k+1);
%%
%% BUILDING KRIGING
%%
%% OBJ 0 BUILDING KGRIGING
% Construct covariance matrix
buildcokriging
%% OBJ 1 BUILDING KGRIGING
% Construct covariance matrix
buildcokriging1
%% OBJ 2 BUILDING KGRIGING
% Construct covariance matrix
buildcokriging2
%%
%%
%
% Make predictions in range 0,1
Xplot=0:0.01:1;
ModelInfo.Option='Pred';
for i=1:101
    pred(i)=cokrigingpredictor(Xplot(i));
    truee(i)=onevar(Xplot(i));
    truec(i)=cheaponevar(Xplot(i));
end

plot(Xplot,truee,'k','LineWidth',2)
hold on
plot(Xplot,truec,'k--','LineWidth',2)
plot(ModelInfo.Xe,ModelInfo.ye,'ko')
plot(ModelInfo.Xc,ModelInfo.yc,'ko')
plot(Xplot,pred,'r')
%legend('f_e','f_c','y_e','y_c','co-Kriging',2)
xlabel('x','FontSize',14)
ylabel('f','FontSize',14)
set(gca,'FontSize',14)
%%
ModelInfo1.Option='Pred';
for i=1:101
    pred1(i)=cokrigingpredictor1(Xplot(i));
    truee(i)=onevar(Xplot(i));
    truec(i)=cheaponevar(Xplot(i));
end

plot(Xplot,truee,'k','LineWidth',2)
hold on
plot(Xplot,truec,'k--','LineWidth',2)
plot(ModelInfo1.Xe,ModelInfo1.ye,'ko')
plot(ModelInfo1.Xc,ModelInfo1.yc,'ko')
plot(Xplot,pred1,'r')
%legend('f_e','f_c','y_e','y_c','co-Kriging',2)
xlabel('x','FontSize',14)
ylabel('f1','FontSize',14)
set(gca,'FontSize',14)
%%
ModelInfo2.Option='Pred';
for i=1:101
    pred2(i)=cokrigingpredictor1(Xplot(i));
    truee(i)=onevar(Xplot(i));
    truec(i)=cheaponevar(Xplot(i));
end
plot(Xplot,truee,'k','LineWidth',2)
hold on
plot(Xplot,truec,'k--','LineWidth',2)
plot(ModelInfo2.Xe,ModelInfo2.ye,'ko')
plot(ModelInfo2.Xc,ModelInfo2.yc,'ko')
plot(Xplot,pred2,'r')
% %legend('f_e','f_c','y_e','y_c','co-Kriging',2)
xlabel('x','FontSize',14)
ylabel('f2','FontSize',14)
set(gca,'FontSize',14)

%%
nV=1; % Number of design variables.
Lb=[0.05]; % Lower bounds of design variables.
Ub=[1]; % Upper bounds of design variables.
% Define parameters of BB-BB algorithm.
nP=100;

OPT=BB_BC(nV,Lb,Ub,nP)




