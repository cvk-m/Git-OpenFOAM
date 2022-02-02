% Calculate the Center of Mass (CM) based on the Big Crunch phase.
function CM=Big_Crunch(P,PFit)
for i=1:size(P,2)
CM(i)=sum(P(:,i)'./PFit)/sum(PFit.^-1);
end