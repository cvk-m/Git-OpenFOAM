function ModelInfo3=buildcokriging3(ModelInfo3)
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
global ModelInfo3
Xe=ModelInfo3.Xe;
Xc=ModelInfo3.Xc;
ye=ModelInfo3.ye;
yc=ModelInfo3.yc;
ne=size(Xe,1); 
nc=size(Xc,1); 
thetad=10.^ModelInfo3.Thetad;
thetac=10.^ModelInfo3.Thetac;
p=1.99;  % added p definition (February 10)
rho=ModelInfo3.rho;

one=ones(ne+nc,1);
y=[yc; ye];
PsicXc=zeros(nc,nc);
for i=1:nc
	for j=i+1:nc
		PsicXc(i,j)=exp(-sum(thetac.*abs(Xc(i,:)-Xc(j,:)).^p));
	end
end
ModelInfo3.PsicXc=PsicXc+PsicXc'+eye(nc)+eye(nc).*eps; 
ModelInfo3.UPsicXc=chol(ModelInfo3.PsicXc);

PsicXe=zeros(ne,ne);
for i=1:ne
	for j=i+1:ne
		PsicXe(i,j)=exp(-sum(thetac.*abs(Xe(i,:)-Xe(j,:)).^p));
	end
end
ModelInfo3.PsicXe=PsicXe+PsicXe'+eye(ne)+eye(ne).*eps; 
ModelInfo3.UPsicXe=chol(ModelInfo3.PsicXe);


PsicXcXe=zeros(nc,ne);
for i=1:nc
	for j=1:ne
		PsicXcXe(i,j)=exp(-sum(thetac.*abs(Xc(i,:)-Xe(j,:)).^p));
	end
end
ModelInfo3.PsicXcXe=PsicXcXe; % Deleted "+[zeros(nc-ne,ne);eye(ne)].*eps" November 2010
ModelInfo3.PsicXeXc=ModelInfo3.PsicXcXe';


PsidXe=zeros(ne,ne);
for i=1:ne
	for j=i+1:ne
		PsidXe(i,j)=exp(-sum(thetad.*abs(Xe(i,:)-Xe(j,:)).^p));
	end
end
ModelInfo3.PsidXe=PsidXe+PsidXe'+eye(ne)+eye(ne).*eps;
ModelInfo3.UPsidXe=chol(ModelInfo3.PsidXe);


ModelInfo3.muc=(ones(nc,1)'*(ModelInfo3.UPsicXc\(ModelInfo3.UPsicXc'\yc)))/(ones(nc,1)'*(ModelInfo3.UPsicXc\(ModelInfo3.UPsicXc'\ones(nc,1))));
ModelInfo3.d=ye-rho.*yc(end-ne+1:end);
ModelInfo3.mud=(ones(ne,1)'*(ModelInfo3.UPsidXe\(ModelInfo3.UPsidXe'\ModelInfo3.d)))/(ones(ne,1)'*(ModelInfo3.UPsidXe\(ModelInfo3.UPsidXe'\ones(ne,1))));

ModelInfo3.SigmaSqrc=(yc-ones(nc,1).*ModelInfo3.muc)'*(ModelInfo3.UPsicXc\(ModelInfo3.UPsicXc'\(yc-ones(nc,1).*ModelInfo3.muc)))/nc; 
ModelInfo3.SigmaSqrd=(ModelInfo3.d-ones(ne,1).*ModelInfo3.mud)'*(ModelInfo3.UPsidXe\(ModelInfo3.UPsidXe'\(ModelInfo3.d-ones(ne,1).*ModelInfo3.mud)))/ne; 

ModelInfo3.C=[ModelInfo3.SigmaSqrc*ModelInfo3.PsicXc rho*ModelInfo3.SigmaSqrc*ModelInfo3.PsicXcXe;
rho*ModelInfo3.SigmaSqrc*ModelInfo3.PsicXeXc rho^2*ModelInfo3.SigmaSqrc*ModelInfo3.PsicXe+ModelInfo3.SigmaSqrd*ModelInfo3.PsidXe];
ModelInfo3.UC=chol(ModelInfo3.C);

ModelInfo3.mu=(one'*(ModelInfo3.UC\(ModelInfo3.UC'\y)))/(one'*(ModelInfo3.UC\(ModelInfo3.UC'\one)));

