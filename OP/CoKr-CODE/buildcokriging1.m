function ModelInfo1=buildcokriging1(ModelInfo1)
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
global ModelInfo1
Xe=ModelInfo1.Xe;
Xc=ModelInfo1.Xc;
ye=ModelInfo1.ye;
yc=ModelInfo1.yc;
ne=size(Xe,1); 
nc=size(Xc,1); 
thetad=10.^ModelInfo1.Thetad;
thetac=10.^ModelInfo1.Thetac;
p=1.99;  % added p definition (February 10)
rho=ModelInfo1.rho;

one=ones(ne+nc,1);
y=[yc; ye];
PsicXc=zeros(nc,nc);
for i=1:nc
	for j=i+1:nc
		PsicXc(i,j)=exp(-sum(thetac.*abs(Xc(i,:)-Xc(j,:)).^p));
	end
end
ModelInfo1.PsicXc=PsicXc+PsicXc'+eye(nc)+eye(nc).*eps; 
ModelInfo1.UPsicXc=chol(ModelInfo1.PsicXc);

PsicXe=zeros(ne,ne);
for i=1:ne
	for j=i+1:ne
		PsicXe(i,j)=exp(-sum(thetac.*abs(Xe(i,:)-Xe(j,:)).^p));
	end
end
ModelInfo1.PsicXe=PsicXe+PsicXe'+eye(ne)+eye(ne).*eps; 
ModelInfo1.UPsicXe=chol(ModelInfo1.PsicXe);


PsicXcXe=zeros(nc,ne);
for i=1:nc
	for j=1:ne
		PsicXcXe(i,j)=exp(-sum(thetac.*abs(Xc(i,:)-Xe(j,:)).^p));
	end
end
ModelInfo1.PsicXcXe=PsicXcXe; % Deleted "+[zeros(nc-ne,ne);eye(ne)].*eps" November 2010
ModelInfo1.PsicXeXc=ModelInfo1.PsicXcXe';


PsidXe=zeros(ne,ne);
for i=1:ne
	for j=i+1:ne
		PsidXe(i,j)=exp(-sum(thetad.*abs(Xe(i,:)-Xe(j,:)).^p));
	end
end
ModelInfo1.PsidXe=PsidXe+PsidXe'+eye(ne)+eye(ne).*eps;
ModelInfo1.UPsidXe=chol(ModelInfo1.PsidXe);


ModelInfo1.muc=(ones(nc,1)'*(ModelInfo1.UPsicXc\(ModelInfo1.UPsicXc'\yc)))/(ones(nc,1)'*(ModelInfo1.UPsicXc\(ModelInfo1.UPsicXc'\ones(nc,1))));
ModelInfo1.d=ye-rho.*yc(end-ne+1:end);
ModelInfo1.mud=(ones(ne,1)'*(ModelInfo1.UPsidXe\(ModelInfo1.UPsidXe'\ModelInfo1.d)))/(ones(ne,1)'*(ModelInfo1.UPsidXe\(ModelInfo1.UPsidXe'\ones(ne,1))));

ModelInfo1.SigmaSqrc=(yc-ones(nc,1).*ModelInfo1.muc)'*(ModelInfo1.UPsicXc\(ModelInfo1.UPsicXc'\(yc-ones(nc,1).*ModelInfo1.muc)))/nc; 
ModelInfo1.SigmaSqrd=(ModelInfo1.d-ones(ne,1).*ModelInfo1.mud)'*(ModelInfo1.UPsidXe\(ModelInfo1.UPsidXe'\(ModelInfo1.d-ones(ne,1).*ModelInfo1.mud)))/ne; 

ModelInfo1.C=[ModelInfo1.SigmaSqrc*ModelInfo1.PsicXc rho*ModelInfo1.SigmaSqrc*ModelInfo1.PsicXcXe;
rho*ModelInfo1.SigmaSqrc*ModelInfo1.PsicXeXc rho^2*ModelInfo1.SigmaSqrc*ModelInfo1.PsicXe+ModelInfo1.SigmaSqrd*ModelInfo1.PsidXe];
ModelInfo1.UC=chol(ModelInfo1.C);

ModelInfo1.mu=(one'*(ModelInfo1.UC\(ModelInfo1.UC'\y)))/(one'*(ModelInfo1.UC\(ModelInfo1.UC'\one)));

