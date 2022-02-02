function OPT=BB_BC(nV,Lb,Ub,nP)
%% Initialization
% Define the properties of COP (tension/compression spring design problem).

% Define parameters of BB-BB algorithm.

maxNFEs=2000; % Maximum number of Objective Function Evaluations.
beta=0.2; % Parameter for controlling the influence of the weighted averageof the particles positions or the center of mass (CM) and the bestparticle.
alfa=1; % Parameter for limiting the size of the initial search space.
%%
%Generate random initial solutions.
%% 
for i=1:nP
P(i,:)=Lb+(Ub-Lb).*rand(1,nV); %Particles matrix or matrix of the initial candidate solutions or the initial population.
end
% Evaluate initial population (P) calling the fobj function constructed inthe second chapter and form its corresponding vectors of objective function(Fit) and penalized objective function (PFit). It should be noted that thedesign vectors all are inside the search space.
for i=1:size(P,1)
[X,fit,pfit]=fobj(P(i,:),Lb,Ub);
P(i,:)=X;
Fit(i)=fit;
PFit(i)=pfit;
end

%%
%%
%Monitor the best candidate solution (bestP) and its correspondingpenalized objective function (minPFit) and objective function (minFit).
%% OBJE 0 BESTpoints
[minPFit,m]=min(PFit);
minFit=Fit(m);
bestP=P(m,:);
%%
%%
%Algorithm Body
NFEs=0; % Current number of Objective Function Evaluations used by thealgorithm until yet.
NITs=0; % Number of algorithm iterations
while NFEs<maxNFEs
% BestP is the best particle found by the BB-BC until yet. MinPFit andMinFit are its corresponding penalized objective function and objectivefunction, respectively.
if NITs==0
BestP=bestP;
MinPFit=minPFit;
MinFit=minFit;
end
NITs=NITs+1; % Update the number of algorithm iterations.
% Generate center of mass vector (CM) based on the Big_Crunch phase.
CM=Big_Crunch(P,PFit);
% Update the particles position and generate new set of solutions basedon the Big_Bang phase.
newP=Big_Bang(P,CM,bestP,beta,alfa,Lb,Ub,NITs);
% Evaluate the new particles. It should be noted that in the BB-BCalgorithm the replacement strategy is not used.
for i=1:size(P,1)
[X,fit,pfit]=fobj(newP(i,:),Lb,Ub);
P(i,:)=X;
Fit(i)=fit;
PFit(i)=pfit;
end
% Update the number of Objective Function Evaluations used by thealgorithm until yet.
NFEs=NFEs+nP;
% Monitor the best candidate solution (bestP) and its correspondingpenalized objective function (minPFit) and objective function (minFit).
[minPFit,m]=min(PFit);
minFit=Fit(m);
bestP=P(m,:);
if minPFit<=MinPFit
BestP=bestP;
MinPFit=minPFit;
MinFit=minFit;
end
% Display desired information of the iteration.
disp(['NITs= ' num2str(NITs) '; MinFit = ' num2str(MinFit) '; MinPFit= ' num2str(MinPFit)]);
% Save the required results for post processing and visualization ofalgorithm performance
output1(NITs,:)=[minFit,minPFit,NFEs];
output2(NITs,:)=[minPFit,max(PFit),mean(PFit)];
output4(NITs,:)=[MinFit,MinPFit,NFEs];
output5(NITs,:)=[BestP,NFEs];
end
%% Monitoring the results
figure(2);
plot((1:1:NITs),output2(:,1),'g',(1:1:NITs),output2(:,2),'r--',(1:1:NITs),output2(:,3),'b-.')
legend('min','max','mean');
xlabel('NITs');
ylabel('pfit');
figure(3);
plot((1:1:NITs),output4(:,2),'g')
xlabel('NITs');
ylabel('MinFIt');
OPT=BestP;