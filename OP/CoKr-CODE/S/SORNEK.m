% Make ModelInfo visible to all functions
global ModelInfo
% Sampling plan
ModelInfo.Xi = Sbestlh(3,12,100,50);
% Compute objective function values – in this case using
% % the dome.m test function
% n – number of points required
% %
% k – number of design variables
% %
% Population – number of individuals in the evolutionary
% %
% operation optimizer
% %
% Iterations – number of generations the evolutionary
% (continued)Sampling Plans
% 27
% %
% operation optimizer is run for
% %
% Note: high values for the two inputs above will ensure high
% %
% quality hypercubes, but the search will take longer.
%ff=[LE,fLl,TP,teta,DA,SW,TE,fT,Rw,wtp,WXA,WYA,AL];

FF=[];
TFP=[];
for i=1:3
FFT=[(0.2*ModelInfo.Xi (i,1)) ,(20*ModelInfo.Xi (i,2)), (ModelInfo.Xi (i,3)),(50*ModelInfo.Xi (i,4))...
    ,(40*ModelInfo.Xi (i,5)-20) ,(40*ModelInfo.Xi (i,6)-20) ,(0.2*ModelInfo.Xi (i,7)),(20*ModelInfo.Xi (i,8))...
    ,(0.15*ModelInfo.Xi (i,9)),(ModelInfo.Xi (i,10) ),(30*ModelInfo.Xi (i,11)-15),(30*ModelInfo.Xi (i,12)-15)]
TFP=cat(1,TFP,FFT);
end 
ModelInfo.X=TFP;












