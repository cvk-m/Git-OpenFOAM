function ModelInfo2=buildcokriging2(ModelInfo2)
% ModelInfo=buildcokriging(ModelInfo)
% Builds co-Kriging model, following parameter estimation. Run prior to
% using cokrigingpredictor.m
%
% Inputs:
%	ModelInfo.Xc - n x k matrix of cheap sample locations
%	ModelInfo.yc - n x 1 vector of cheap observed data
%	ModelInfo.Xe - n x k matrix of expensive sample locations
%	ModelInfo.ye - n x 1 vector of expensive observed data
%   ModelInfo.Thetac - 1 x k vector of log(theta) of cheap model
%   ModelInfo.Thetad - 1 x k vector of log(theta) of difference model
%   ModelInfo.rho - scalar scaling parameter
% Outputs:
%	ModelInfo
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

% global ModelInfo %switch to 'anonymous functions' June 2016
global ModelInfo2
Xe=ModelInfo2.Xe;
Xc=ModelInfo2.Xc;
ye=ModelInfo2.ye;
yc=ModelInfo2.yc;
ne=size(Xe,1); 
nc=size(Xc,1); 
thetad=10.^ModelInfo2.Thetad;
thetac=10.^ModelInfo2.Thetac;
p=1.99;  % added p definition (February 10)
rho=ModelInfo2.rho;

one=ones(ne+nc,1);
y=[yc; ye];
PsicXc=zeros(nc,nc);
for i=1:nc
	for j=i+1:nc
		PsicXc(i,j)=exp(-sum(thetac.*abs(Xc(i,:)-Xc(j,:)).^p));
	end
end
ModelInfo2.PsicXc=PsicXc+PsicXc'+eye(nc)+eye(nc).*eps; 
ModelInfo2.UPsicXc=chol(ModelInfo2.PsicXc);

PsicXe=zeros(ne,ne);
for i=1:ne
	for j=i+1:ne
		PsicXe(i,j)=exp(-sum(thetac.*abs(Xe(i,:)-Xe(j,:)).^p));
	end
end
ModelInfo2.PsicXe=PsicXe+PsicXe'+eye(ne)+eye(ne).*eps; 
ModelInfo2.UPsicXe=chol(ModelInfo2.PsicXe);


PsicXcXe=zeros(nc,ne);
for i=1:nc
	for j=1:ne
		PsicXcXe(i,j)=exp(-sum(thetac.*abs(Xc(i,:)-Xe(j,:)).^p));
	end
end
ModelInfo2.PsicXcXe=PsicXcXe; % Deleted "+[zeros(nc-ne,ne);eye(ne)].*eps" November 2010
ModelInfo2.PsicXeXc=ModelInfo2.PsicXcXe';


PsidXe=zeros(ne,ne);
for i=1:ne
	for j=i+1:ne
		PsidXe(i,j)=exp(-sum(thetad.*abs(Xe(i,:)-Xe(j,:)).^p));
	end
end
ModelInfo2.PsidXe=PsidXe+PsidXe'+eye(ne)+eye(ne).*eps;
ModelInfo2.UPsidXe=chol(ModelInfo2.PsidXe);


ModelInfo2.muc=(ones(nc,1)'*(ModelInfo2.UPsicXc\(ModelInfo2.UPsicXc'\yc)))/(ones(nc,1)'*(ModelInfo2.UPsicXc\(ModelInfo2.UPsicXc'\ones(nc,1))));
ModelInfo2.d=ye-rho.*yc(end-ne+1:end);
ModelInfo2.mud=(ones(ne,1)'*(ModelInfo2.UPsidXe\(ModelInfo2.UPsidXe'\ModelInfo2.d)))/(ones(ne,1)'*(ModelInfo2.UPsidXe\(ModelInfo2.UPsidXe'\ones(ne,1))));

ModelInfo2.SigmaSqrc=(yc-ones(nc,1).*ModelInfo2.muc)'*(ModelInfo2.UPsicXc\(ModelInfo2.UPsicXc'\(yc-ones(nc,1).*ModelInfo2.muc)))/nc; 
ModelInfo2.SigmaSqrd=(ModelInfo2.d-ones(ne,1).*ModelInfo2.mud)'*(ModelInfo2.UPsidXe\(ModelInfo2.UPsidXe'\(ModelInfo2.d-ones(ne,1).*ModelInfo2.mud)))/ne; 

ModelInfo2.C=[ModelInfo2.SigmaSqrc*ModelInfo2.PsicXc rho*ModelInfo2.SigmaSqrc*ModelInfo2.PsicXcXe;
rho*ModelInfo2.SigmaSqrc*ModelInfo2.PsicXeXc rho^2*ModelInfo2.SigmaSqrc*ModelInfo2.PsicXe+ModelInfo2.SigmaSqrd*ModelInfo2.PsidXe];
ModelInfo2.UC=chol(ModelInfo2.C);

ModelInfo2.mu=(one'*(ModelInfo2.UC\(ModelInfo2.UC'\y)))/(one'*(ModelInfo2.UC\(ModelInfo2.UC'\one)));

