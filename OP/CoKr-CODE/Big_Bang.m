% Generate new set of particles based on the Big_Bang phase
function newP=Big_Bang(P,CM,bestP,beta,alfa,Lb,Ub,NITs)
for i=1:size(P,1)
newP(i,:)=beta*CM+(1-beta)*bestP+randn(1,size(P,2)).*((Ub-Lb)/(NITs));
end