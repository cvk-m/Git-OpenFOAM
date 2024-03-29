clc
clear all
DM=0; % if 3d mesh =1 ,0 means surface mesh
Nfoil=200;% number of point to built airfoil base shape
Nf=Nfoil+1;
dbeta=1/Nfoil;
LE=0;  %LEADING EDGE MODIFICATION  LENGHT  
fLl = 0; %frequency
Nhp=300; % 3d için katman sayısı
TP=1;
AL=1;
teta =0; %randi([2 50],1,1)twist
DA=0;
SW=0;
Nw=0;
TE=0;  %TRAILING EDGE MODIFICATION  LENGHT  
fT=10;
Rw=0.125;
wtp=1;  %taper
WXA=10;
WYA=10;
fL = fT;
ff=[LE,fLl,TP,teta,DA,SW,TE,fT,Rw,wtp,WXA,WYA,AL];
fWe=[];
for i=1:13
    if i==13

fcx=mat2str(ff(i));
fWe=cat(2,fWe,fcx); 
    else
 fxx='-';
fx=mat2str(ff(i));
fcx=strcat(fx,fxx);
fWe=cat(2,fWe,fcx);
    end
end
f11='.geo';

f13='w';

fWE=strcat(fWe,f11);
fid = fopen(fWE,f13);
%fid = fopen('1.geo','w');


nx=round(Nhp/fT);

while 1
    Nax=fT*nx;
if Nax>= Nhp
     nm=mod(Nax,2);
if nm==0
    Na=Nax+1;
else
Na=fT*(nx+1)+1; % NUMBER OF AIRFOL TO BUILT 3D WING    
end
break
end
nx=nx+1;
end 


%Na=101;
Ns=Na+Nw;%round((AL*Rw)*Na/(AL));

h=0.1; % mesh size
lms=0.5;
tms=0.5;
k=1;
  % AIRFOIL LENGHT
%%
beta=0;
x = zeros(1,Nf);

for i=0:dbeta:1
   x(k)=i;
    beta=beta+dbeta;
    k=k+1;
end


NACA12= 'NACA0012-201.xlsx';
NACA0012=xlsread(NACA12);
x=NACA0012(1,:);
yufoil=NACA0012(2,:);
ylfoil=NACA0012(3,:);


%% AIRFOIL LENGHT COORD.
a=1;
dNZ=AL/(Na-1);
z = zeros(1,Na);
for i=0:dNZ:AL  %% cox
z(a)=i;
a=a+1;
end
YU=zeros(1,Na*Nf);
YL=zeros(1,Na*Nf);
for i=1:Na
    for a=1:Nf
    YU(a+(i-1)*Nf)=yufoil(a);
    YL(a+(i-1)*Nf)=ylfoil(a);
    end
end
Z=zeros(1,Na*Nf);
for i=1:Na
    for a=(i-1)*Nf+1:1:(Nf*i)
    Z(a)=z(i);
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LEADING EDGE EXTENTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BASE SIN WAVE FUNC FOR LEADING


LEP=[];
zL=[];

for aa=0:dNZ:AL
 le=LE/2 ;
zQL =le*cos(pi * fLl * aa+pi)-le;
LEP=cat(1,LEP,-zQL);
zL=cat(1,zL,aa);
end

%plot(zi,LEP)
%% LEADING EDGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFICATION
[idxu]=max(yufoil(:));
Yu=min(find(yufoil==idxu));
%[max_num, max_idx]=max(yufoil(:))
%[X,Y]=ind2sub(size(yufoil),find(yufoil==max_num))  % Y number of point to adjust new tubcle shape airfoil
Lxu=[];
s=1;
LXu=[];
ru=zeros(1,Yu-1);

for c=1:Na%NZ-1   %UPPER SURFACE MODIFIED
stb=LEP(c); %extent max. span of tubcle

if stb==0
   ru(1:Yu-1)=x(1:Yu-1);
else
dbu=(stb+x(Yu-1))/(Yu-2);
beta=0;
for g=1:(Yu-1)
   ru(g)=-stb+beta;
    beta=beta+dbu;
end
end
ru(Yu-1)=(ru(Yu-1)+x(Yu-1))/2;
Lxu=[ru,x(1,(Yu):Nf)];
LXu=cat(2,LXu,Lxu);
end
[idxl]=min(ylfoil(:)) ;   %LOWER SURFACE MODIFIED
Yl=min(find(ylfoil==idxl));
Lxl=[];
LXl=[];
rl=zeros(1,Yl-1);


for c=1:Na %NZ-1
  stb=LEP(c) ;
beta=0;
if stb==0

  rl(1:Yl-1)=x(1:Yl-1);

else
dbl=(stb+x(Yl-1))/(Yl-2);
beta=0;
for g=1:(Yl-1)
   rl(g)=-stb+beta;
    beta=beta+dbl;
end
end
rl(Yl-1)=(rl(Yl-1)+x(Yl-1))/2;
Lxl=[rl,x(1,(Yl):Nf)];
LXl=cat(2,LXl,Lxl);
end
%LXl=-(((-LXl(1,1)+LXl(1,Nf))/2)+LXl(1,1))+LXl;
%LXu=-(((-LXu(1,1)+LXu(1,Nf))/2)+LXu(1,1))+LXu;

%% TRAILING SHARP EDGE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BASE SIN WAVE FUNC FOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRAILING EDGE

TEP=[];
zT=[];

for aa=0:dNZ:AL
 tek=TE/2  ;
zQT =tek*cos(pi * fT * aa)-tek;
TEP=cat(1,TEP,-zQT);
zT=cat(1,zT,aa);
end

%%
  TXu=LXu;
  NN=[];
  
 for i=1:Na %NZ-1
      b=1;
     for a=(i-1)*Nf+1:1:(Nf*i) 
if TEP(i)==0
    TXu(a)=TXu(a);
    NN(i,b)=0;
  b=b+1;
elseif LXu(a)>((0.5+LE/2)-TEP(i))
   TXu(a)=((0.5+LE/2)-TEP(i));
NN(i,b)=a;
  b=b+1;
else
    TXu(a)=TXu(a);
end
     end
 end

    TXl=TXu;

%end

uTP=zeros(1,Na*Nf);
xTP=zeros(1,Na*Nf);
lTP=zeros(1,Na*Nf)  ;%y duzenlenecek
xTPl=zeros(1,Na*Nf);
for i=1:Na %NZ
   dTp=(1-TP)/(Na-1);
    Tp=1-dTp*(i-1);
    for a=(i-1)*Nf+1:1:(Nf*i)
      uTP(a)=YU(a)*Tp  ;%y duzenlenecek
      xTP(a)=TXu(a)*Tp ; %x duzenlenecek
      lTP(a)=YL(a)*Tp  ;%y duzenlenecek
      xTPl(a)=TXl(a)*Tp;
     end
 end
%%   twst  
UXTP=transpose([xTP;uTP;Z]);
LXTP=transpose([xTPl;lTP;Z]);
%PU=transpose([TXu+0.1;YU;Z])

sw=SW*pi/180;
te=teta*pi/180;
da=DA*pi/180;
LTWW=[];
TWW=[];
for i=1:Na
    t=te/(Na-1);
    AT=t*(i-1);
    r=rotz(AT);
    TW=UXTP((i-1)*Nf+1:(i*Nf),:)*r;
    LTW=LXTP((i-1)*Nf+1:(i*Nf),:)*r;
    TWW=cat(1,TWW,TW);
    LTWW=cat(1,LTWW,LTW);
end
%%  TAPER
%WPU=[]
%for i=1:2
    %Tp=1+((TP-1)*(i-1)/2)
    %WP=[TWW((i-1)*201+1:(i*201),1:2),PU((i-1)*201+1:(i*201),3)]
    %WPU=cat(1,WPU,WP);
%end
   
WUXTP=transpose([xTP;uTP;Z]);
%% SWEEP
LSWPU=[];
SWPU=[];
for i=1:Na
    SS=sw*(i-1)/(Na-1);
    SWP=[TWW((i-1)*Nf+1:(i*Nf),1)+transpose(Z((i-1)*Nf+1:(i*Nf))*tan(SS)),TWW((i-1)*Nf+1:(i*Nf),2),transpose(Z((i-1)*Nf+1:(i*Nf)))];
    SWPU=cat(1,SWPU,SWP);
    LSWP=[LTWW((i-1)*Nf+1:(i*Nf),1)+transpose(Z((i-1)*Nf+1:(i*Nf))*tan(SS)),LTWW((i-1)*Nf+1:(i*Nf),2),transpose(Z((i-1)*Nf+1:(i*Nf)))];
    LSWPU=cat(1,LSWPU,LSWP);
end
   
%%  
  LDSWPU=[];
  DSWPU=[];
for i=1:Na
     DS=da*(i-1)/(Na-1);
     DWP=[SWPU((i-1)*Nf+1:(i*Nf),1),SWPU((i-1)*Nf+1:(i*Nf),2)-transpose(Z((i-1)*Nf+1:(i*Nf))*tan(DS)),transpose(Z((i-1)*Nf+1:(i*Nf)))];
     DSWPU=cat(1,DSWPU,DWP);
     LDWP=[LSWPU((i-1)*Nf+1:(i*Nf),1),LSWPU((i-1)*Nf+1:(i*Nf),2)-transpose(Z((i-1)*Nf+1:(i*Nf))*tan(DS)),transpose(Z((i-1)*Nf+1:(i*Nf)))];
     LDSWPU=cat(1,LDSWPU,LDWP);
end




%% WINGLET


wax=WXA*pi/180;
way=WYA*pi/180;
WL=AL*Rw/Nw;
wwtp=(1-wtp)/(Nw-1);
WTA=[];
 LWTA=[];
 if Nw==0
     ARFLWu=DSWPU;
     ARFLWl=LDSWPU;
     ARFLW=[DSWPU;LDSWPU];
 else
for i=1:Nw
    WP=[DSWPU((Na-1)*Nf+1:Na*Nf,1),DSWPU((Na-1)*Nf+1:Na*Nf,2),DSWPU((Na-1)*Nf+1:Na*Nf,3)+WL*i];
    WA=[WP(:,1)-WL*i*tan(wax),WP(:,2)-WL*i*tan(way),WP(:,3)];
    WAT=[WA(:,1)*(1-i*wwtp),WA(:,2)*(1-i*wwtp),WA(:,3)];
    WTA=cat(1,WTA,WAT);
   
    LWP=[LDSWPU((Na-1)*Nf+1:Na*Nf,1),LDSWPU((Na-1)*Nf+1:Na*Nf,2),LDSWPU((Na-1)*Nf+1:Na*Nf,3)+WL*i];
    LWA=[LWP(:,1)-WL*i*tan(wax),LWP(:,2)-WL*i*tan(way),LWP(:,3)];
    LWAT=[LWA(:,1)*(1-i*wwtp),LWA(:,2)*(1-i*wwtp),LWA(:,3)];
    LWTA=cat(1,LWTA,LWAT);
   
end
ARFLWu=[DSWPU;WTA];
ARFLWl=[LDSWPU;LWTA];
ARFLW=[ARFLWu;ARFLWl];
 end
      ARFLWuy=[];
     ARFLWly=[];
 inx=1;
for d=1:Ns
 for i=1:Nf-1
   ARFLWuy(inx,1) = ARFLWu(i+(d-1)*Nf,1);
   ARFLWuy(inx,2) = ARFLWu(i+(d-1)*Nf,2);
   ARFLWuy(inx,3) = ARFLWu(i+(d-1)*Nf,3);
   
   ARFLWly(inx,1) = ARFLWl(i+(d-1)*Nf+1,1);
   ARFLWly(inx,2) = ARFLWl(i+(d-1)*Nf+1,2);  
   ARFLWly(inx,3) = ARFLWl(i+(d-1)*Nf+1,3);
   inx=inx+1;
 end
end
ARFLWy=[ARFLWuy;ARFLWly];
%%%%% FİLE NAME mesh size 
    [numRows,numCols]=size(ARFLWy);
ff=h*ones(numRows,1);
ff1=ff;
for i=1:Ns
    for a=1:Nf
 if a== Nf
        ff1(a+(i-1)*2*(Nf-1))=ff(a+(i-1)*2*(Nf-1))*tms;
 elseif a==1
         ff1(a+(i-1)*2*(Nf-1))=ff(a+(i-1)*2*(Nf-1))*lms;
        elseif a <= (Nf-1)*0.36
        ff1(a+(i-1)*2*(Nf-1))=ff(a+(i-1)*2*(Nf-1))*lms;
        ff1(2*(Nf+(2*(i-1)*(Nf-1)))-(a+(i-1)*2*(Nf-1)))= ff(2*(Nf+(2*(i-1)*(Nf-1)))-(a+(i-1)*2*(Nf-1)))*lms;

        elseif (a>= Nf*0.75)
          ff1(a+(i-1)*2*(Nf-1))=ff(a+(i-1)*2*(Nf-1))*tms; 
          ff1(2*(Nf+(2*(i-1)*(Nf-1)))-(a+(i-1)*2*(Nf-1)))= ff(2*(Nf+(2*(i-1)*(Nf-1)))-(a+(i-1)*2*(Nf-1)))*tms;
        else
            
 end
    end 
end

%%%%% mesh point 
   [r,c]=size(ARFLWy(:,1));
   gridPts=r;
    point = zeros(r,4);
    ind = 1;
    for d = 1:1:Ns
        
    for i = 1:1:Nf-1
        point(ind,1) = ARFLWuy(i+(d-1)*(Nf-1),1);
        point(Nf-1+ind,1) = ARFLWly(Nf+(d-1)*(Nf-1)-i,1);
        
        point(ind,2) = ARFLWuy(i+(d-1)*(Nf-1),2);
        point(Nf-1+ind,2) = ARFLWly(Nf+(d-1)*(Nf-1)-i,2);
        
        point(ind,3) = ARFLWuy(i+(d-1)*(Nf-1),3);
        point(Nf-1+ind,3) = ARFLWly(Nf+(d-1)*(Nf-1)-i,3);
        
        ind = ind + 1;
    end
    ind = ind+Nf-1;
    end
   
   % point(:,4) =h*ones(size(r));
       point(:,4) = ff1;
     assignin('base','point',point);
    
    % Number of points and lines
    numPts = length(point(:,1));
    numLns = numPts;
    
    %% AIRFOIL POINTS

    for i = 1:1:numPts
        fprintf(fid,'Point(%i) = {%g, %g, %g, %g};\r\n',i,point(i,1),point(i,2),...
                                                    point(i,3),point(i,4));
    end

   %%%%%%%%%%%%%%%% surface airfoil front side

taf=Nf-(round((Nf)*0.3)+1);  %
Ntaf=(round((Nf)*0.3)+1);   
 ptafu=[];
 ptafl=[];
for i=1:Na
    ptafu(i)=Nf+(i-1)*2*(Nf-1)-Ntaf;
    ptafl(i)=Nf+(i-1)*2*(Nf-1)+Ntaf;  
end   

%%
%ltaf 

Ltaf=fL;
laF=round(0.3*Nf); %01 idi 03
 pLtafu=[];
 pLtafl=[];
for i=1:Na
pLtafu(i)=(i-1)*(Nf-1)*2+laF;
pLtafl(i)=i*(Nf-1)*2-laF;  
end

%%

x1='%i, ';
x11='%i';
y1='Line(%i) = {';
z1='}';
zx1=';\r\n';
y12='Transfinite Curve {';
z12='}= 200 Using Progression 1}';
%
We=[];

for i=1:Na
if (i==Na)
        k1=x11;
    else 
   k1=x1;
end 
 We=cat(2,We,k1);
end
a=0;    
a1=[a:(Nf-1)*2:(Nf-1)*2*(Na-1)+1];
spaceFormat1=[y1 We z1 zx1];


%%
% Trailing edge lines
nn=[];
for i=1:Na-1
if NN(i,1)==0
nnxx(i)=Nf+(i-1)*(Nf-1)*2;
nnyy(i)=Nf+(i-1)*(Nf-1)*2;
nn(i)=0;
    else
nn(i)=(i)*Nf-NN(i,1);
nnxx(i)=1+((Nf-1)*(2*i-1))-nn(i); 
nnyy(i)=1+((Nf-1)*(2*i-1))+nn(i);
end
end
nnxx=[nnxx,1+((Nf-1)*(2*Na-1))];
nnyy=[nnyy,1+((Nf-1)*(2*Na-1))];
%% surf lines
for i=1:Na

   if (i==1)
fprintf(fid,'Line(%i) = {%i:%i};\r\n',a+6*(i-1)+1,1, pLtafu(i));
fprintf(fid,'Line(%i) = {%i:%i};\r\n',a+6*(i-1)+2,pLtafu(i),ptafu(i));
fprintf(fid,'Line(%i) = {%i:%i};\r\n',a+6*(i-1)+3,ptafu(i),nnxx(i));
fprintf(fid,'Line(%i) = {%i:%i};\r\n',a+6*(i-1)+4,nnyy(i),ptafl(i));
fprintf(fid,'Line(%i) = {%i:%i};\r\n',a+6*(i-1)+5,ptafl(i),pLtafl(i));
fprintf(fid,'Line(%i) = {%i:%i,%i};\r\n',a+6*(i-1)+6,pLtafl(i),2*(Nf-1),1);
    else
fprintf(fid,'Line(%i) = {%i:%i};\r\n',a+6*(i-1)+1,1+(i-1)*2*(Nf-1),pLtafu(i));
fprintf(fid,'Line(%i) = {%i:%i};\r\n',a+6*(i-1)+2,pLtafu(i),ptafu(i));
fprintf(fid,'Line(%i) = {%i:%i};\r\n',a+6*(i-1)+3,ptafu(i),nnxx(i));%Nf+(i-1)*2*(Nf-1)
fprintf(fid,'Line(%i) = {%i:%i};\r\n',a+6*(i-1)+4,nnyy(i),ptafl(i));%Nf+(i-1)*2*(Nf-1)
fprintf(fid,'Line(%i) = {%i:%i};\r\n',a+6*(i-1)+5,ptafl(i),pLtafl(i));%Nf+(i-1)*2*(Nf-1)
fprintf(fid,'Line(%i) = {%i:%i,%i};\r\n',a+6*(i-1)+6,pLtafl(i),((i-1)*2+2)*(Nf-1),(i-1)*2*(Nf-1)+1);
    end
end


Ww=[];
for i=1:(Na-1)/fT+1
if (i==(Na-1)/fT+1)
        k1=x11; %burada deneme x11 düzelt unutursan
    else 
   k1=x1;
end 
 Ww=cat(2,Ww,k1);
end

% kesik üçgen iki parça + bölen çizgi iptal
spaceFormat2=[y1 Ww z1 zx1];

% kanat sonu her slot için kenar çizgisi 
 Ll=round((Na-1)/fL)-1;
  L2=round((Na-1)/fL)+1;
  L0=round((Na-1)/fL);
 W4=[];
for i=1:L0
if (i==L0)
        k1=x11; %burada deneme x11 düzelt unutursan
    else 
   k1=x1;
end 
 W4=cat(2,W4,k1);
end
W5=[];
for i=(fL-1)*L0+1:1:Na
if (i==Na)
        k1=x11; %burada deneme x11 düzelt unutursan
    else 
   k1=x1;
end 
 W5=cat(2,W5,k1);
end
W6=[];
for i=1:L2
if (i==L2)
        k1=x11; %burada deneme x11 düzelt unutursan
    else 
   k1=x1;
end 
 W6=cat(2,W6,k1);
end
spaceFormat4=[y1 W4 z1 zx1];
spaceFormat5=[y1 W5 z1 zx1];
spaceFormat6=[y1 W6 z1 zx1];

s1=a+6*(Na-1)+6;

% tra1   s1+i s1+fT+i  
%%
for i=1:fL
    if i==fL
fprintf(fid,spaceFormat5,s1+i,nnxx((i-1)*(L0)+1:Na));   
    elseif i==1
fprintf(fid,spaceFormat6,s1+i,nnxx(1:(L0+1)));        
    else   
fprintf(fid,spaceFormat6,s1+i,nnxx((i-1)*(L0)+1:i*(L0)+1));
    end
end

for i=1:fL
     if i==fL
         
fprintf(fid,spaceFormat5,s1+fT+i,nnyy((i-1)*(L0)+1:Na));     
    elseif i==1
 fprintf(fid,spaceFormat6,s1+fT+i,nnyy(1:(L0+1)));    
     
    else 
fprintf(fid,spaceFormat6,s1+fT+i,nnyy((i-1)*(L0)+1:i*(L0)+1));
      end
end

%+ bölen çizgi b
s2=s1+2*fT;
if round(fT/2)-fT/2==0
for i=1:round(fT/2)
if i==1
fprintf(fid,'Line(%i) = {%i,%i};\r\n',s2+i,nnxx(L0+1),nnyy(L0+1));
else
fprintf(fid,'Line(%i) = {%i,%i};\r\n',s2+i,nnxx(L0*(2*i-1)+1),nnyy(L0*(2*i-1)+1));
end
end
else
for i=1:round(fT/2)
if i==1
fprintf(fid,'Line(%i) = {%i,%i};\r\n',s2+i,nnxx(L0+1),nnyy(L0+1));
elseif i==round(fT/2)
fprintf(fid,'Line(%i) = {%i,%i};\r\n',s2+i,nnxx(Na),nnyy(Na));
end
end
end
s3=s2+round(fT/2);
% ntaf 

for i=1:fL
     if i==fL
         
fprintf(fid,spaceFormat5,s3+i,ptafu((i-1)*L0+1:Na));     
    elseif i==1
 fprintf(fid,spaceFormat6,s3+i,ptafu(1:i*L0+1));    
     
    else 
fprintf(fid,spaceFormat6,s3+i,ptafu((i-1)*L0+1:i*L0+1));
      end
end

s4=s3+fT;
for i=1:fL
     if i==fL
         
fprintf(fid,spaceFormat5,s4+i,ptafl((i-1)*L0+1:Na));     
    elseif i==1
 fprintf(fid,spaceFormat6,s4+i,ptafl(1:i*L0+1));    
     
    else 
fprintf(fid,spaceFormat6,s4+i,ptafl((i-1)*L0+1:i*L0+1));
      end
end
%ntaf ile son çizgi arası üst yüzey bölme nokta+ son çizgi ayıran nokta  k1
%k2


s5=s4+fT;

%WİNGLET LİNE
%leading edge line
W3=[];
s7=s5;

%fprintf(fid,'Line(%i) = {%i: %i};\r\n',s7+1,1,Nf);
%fprintf(fid,'Line(%i) = {%i: %i, %i};\r\n',s7+2,Nf,2*(Nf-1),1);
s8=s7+2;
if Nw==0
fprintf(fid,'Line(%i) = {%i: %i};\r\n',s8+1,6*(Na-1)+1,6*(Na-1)+Nf);
fprintf(fid,'Line(%i) = {%i: %i, %i};\r\n',s8+2,6*(Na-1)+Nf,6*(Na-1)+2*(Nf-1),6*(Na-1)+1);
s9=s8+2;
else
for i=Na:Ns
if (i==Ns)
        k1=x11; %burada deneme x11 düzelt unutursan
    else 
   k1=x1;
end 
 W3=cat(2,W3,k1);
end
spaceFormat3=[y1 W3 z1 zx1];
aw3=((Na-1)*2*(Nf-1)+1:(Nf-1)*2:(Ns-1)*2*(Nf-1)+1);
fprintf(fid,spaceFormat3,s8+1,aw3);
%trailing edge
aw3t=((Na-1)*2*(Nf-1)+1+(Nf-1):(Nf-1)*2:(Ns-1)*2*(Nf-1)+Nf);
fprintf(fid,spaceFormat3,s8+2,aw3t);
%firts winglet airfoil edge up

fprintf(fid,'Line(%i) = {%i: %i};\r\n',s8+3,(Ns-1)*2*(Nf-1)+1,(Ns-1)*2*(Nf-1)+laF);
fprintf(fid,'Line(%i) = {%i: %i};\r\n',s8+4,(Ns-1)*2*(Nf-1)+laF,(Ns-1)*2*(Nf-1)+Nf-Ntaf);
fprintf(fid,'Line(%i) = {%i: %i};\r\n',s8+5,(Ns-1)*2*(Nf-1)+Nf-Ntaf,(Ns-1)*2*(Nf-1)+Nf);
wpLu=(Na-1)*2*(Nf-1)+laF:2*(Nf-1):(Ns-1)*2*(Nf-1)+laF;
wpNu=(Na-1)*2*(Nf-1)+Nf-Ntaf:2*(Nf-1):(Ns-1)*2*(Nf-1)+Nf-Ntaf;
fprintf(fid,spaceFormat3,s8+6,wpLu);
fprintf(fid,spaceFormat3,s8+7,wpNu);
%fprintf(fid,'Line(%i) = {%i: %i};\r\n',s8+1,(Ns-1)*2*(Nf-1)+1,(Ns-1)*2*(Nf-1)+Nf);
%last winglet airfoil edge low
fprintf(fid,'Line(%i) = {%i: %i};\r\n',s8+8,(Ns-1)*2*(Nf-1)+Nf,(Ns-1)*2*(Nf-1)+Nf+Ntaf);
fprintf(fid,'Line(%i) = {%i: %i};\r\n',s8+9,(Ns-1)*2*(Nf-1)+Nf+Ntaf,(Ns-1)*2*(Nf-1)+2*(Nf-1)-laF);
fprintf(fid,'Line(%i) = {%i: %i, %i};\r\n',s8+10,(Ns-1)*2*(Nf-1)+2*(Nf-1)-laF,(Ns-1)*2*(Nf-1)+2*(Nf-1),(Ns-1)*2*(Nf-1)+1);
%fprintf(fid,'Line(%i) = {%i: %i, %i};\r\n',s8+2,(Ns-1)*2*(Nf-1)+Nf,(Ns-1)*2*(Nf-1)+2*(Nf-1),(Ns-1)*2*(Nf-1)+1);

wpLl=(Na-1)*2*(Nf-1)+2*(Nf-1)-laF:2*(Nf-1):(Ns-1)*2*(Nf-1)+2*(Nf-1)-laF;
wpNl=(Na-1)*2*(Nf-1)+Nf+Ntaf:2*(Nf-1):(Ns-1)*2*(Nf-1)+Nf+Ntaf;

fprintf(fid,spaceFormat3,s8+11,wpNl);
fprintf(fid,spaceFormat3,s8+12,wpLl);
%fprintf(fid,'Line(%i) = {%i: %i};\r\n',s8+13,(Ns-1)*2*(Nf-1)+1,(Ns-1)*2*(Nf-1)+Nf);
%fprintf(fid,'Line(%i) = {%i: %i, %i};\r\n',s8+14,(Ns-1)*2*(Nf-1)+Nf,(Ns-1)*2*(Nf-1)+2*(Nf-1),(Ns-1)*2*(Nf-1)+1);


s9=s8+12;
%fprintf(fid,'Line(%i) = {%i: %i};\r\n',s8+1,(Na-1)*2*(Nf-1)+1,(Na-1)*2*(Nf-1)+1+(Nf-1));
end
%firts winglet airfoil edge low
%fprintf(fid,'Line(%i) = {%i: %i, %i};\r\n',s8+2,(Na-1)*2*(Nf-1)+1+(Nf-1),(Na-1)*2*(Nf-1)+2*(Nf-1),(Na-1)*2*(Nf-1)+1);

%last winglet airfoiil edge up

%1 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ÜST YÜZEY ALANI BÖLgE İÇİN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ÇİZGİLER
% leading edge çizgi çizgi
 a1=1:(Nf-1)*2:(Nf-1)*2*(Na-1)+1;



for i=1:fL
    if i==fL
fprintf(fid,spaceFormat5,s9+i,a1((i-1)*(L0)+1:Na));   
    elseif i==1
fprintf(fid,spaceFormat6,s9+i,a1(1:(L0)+1));        
    else   
fprintf(fid,spaceFormat6,s9+i,a1((i-1)*(L0)+1:i*(L0)+1));
    end
end
%s10=s9+fL; %
%for i=1:fL
 %  if i==fL
%fprintf(fid,spaceFormat5,s10+i,ptafu((i-1)*L0:Na));     
   %elseif i==1
   %  fprintf(fid,spaceFormat4,s10+i,ptafu(1:i*L0));
  % else 
  %  fprintf(fid,spaceFormat6,s10+i,ptafu((i-1)*L0:i*L0));
 %  end
%end
%s12=s11+fL;

%for i=1:fL
  %   if i==fL
         
%fprintf(fid,spaceFormat5,s12+i,ptafl((i-1)*L0:Na));     
 %   elseif i==1
% fprintf(fid,spaceFormat4,s12+i,ptafl(1:i*L0));    
     
%    else 
%fprintf(fid,spaceFormat6,s12+i,ptafl((i-1)*L0:i*L0));
  %    end
%end
%%
%ltaf
s10=s9+fL;
for i=1:fL
     if i==fL
         
fprintf(fid,spaceFormat5,s10+i,pLtafu((i-1)*L0+1:Na));     
    elseif i==1
 fprintf(fid,spaceFormat6,s10+i,pLtafu(1:i*L0+1));    
     
    else 
fprintf(fid,spaceFormat6,s10+i,pLtafu((i-1)*L0+1:i*L0+1));
      end
end

s11=s10+fL;
for i=1:fL
     if i==fL
         
fprintf(fid,spaceFormat5,s11+i,pLtafl((i-1)*L0+1:Na));     
    elseif i==1
 fprintf(fid,spaceFormat6,s11+i,pLtafl(1:i*L0+1));    
     
    else 
fprintf(fid,spaceFormat6,s11+i,pLtafl((i-1)*L0+1:i*L0+1));
      end
end
%s12=s11+fL;
%winglet single airfoil line  up
%fprintf(fid,'Line(%i) = {%i: %i};\r\n',s12+1,(Na-1)*2*(Nf-1)+1,(Na-1)*2*(Nf-1)+1+(Nf-1));
%winglet single airfoil line  low
%fprintf(fid,'Line(%i) = {%i: %i, %i};\r\n',s12+2,(Na-1)*2*(Nf-1)+1+(Nf-1),(Na-1)*2*(Nf-1)+2*(Nf-1),(Na-1)*2*(Nf-1)+1);
%%

if (round(fT/2)-fT/2)==0 %%% slot içi üzey 
for b=1:round(fT/2)
fprintf(fid,'Line Loop(%i) = {%i,%i,%i};\r\n',b,s1+2*b-1,-(s1+2*b-1+fT),s1+2*fT+b);
fprintf(fid,'Plane Surface(%i) = {%i};\r\n',b,b);
fprintf(fid,'Line Loop(%i) = {%i,%i,%i};\r\n',round(fT/2)+b,s1+2*b,-(s1+2*b+fT),-(s1+2*fT+b));
fprintf(fid,'Plane Surface(%i) = {%i};\r\n',round(fT/2)+b,round(fT/2)+b);    
end
else
for b=1:round(fT/2)
if b==(round(fT/2))
fprintf(fid,'Line Loop(%i) = {%i,%i,%i};\r\n',b,s1+b,-(s1+b+fT),s1+2*fT+b);
fprintf(fid,'Plane Surface(%i) = {%i};\r\n',b,b);    
else
fprintf(fid,'Line Loop(%i) = {%i,%i,%i};\r\n',b,s1+2*b-1,-(s1+2*b-1+fT),s1+2*fT+b);
fprintf(fid,'Plane Surface(%i) = {%i};\r\n',b,b);
fprintf(fid,'Line Loop(%i) = {%i,%i,%i};\r\n',round(fT/2)+b,s1+2*b,-(s1+2*b+fT),-(s1+2*fT+b));
fprintf(fid,'Plane Surface(%i) = {%i};\r\n',round(fT/2)+b,round(fT/2)+b);
end
end
end

y1=round(fT/2)+b;

if DM==1

% yan yüzeyler
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i,%i,%i};\r\n',y1+1,1,2,3,4,5,6);
%fprintf(fid,'Line Loop(%i) = {%i,%i};\r\n',y1+1,s7+1,s7+2);
fprintf(fid,'Plane Surface(%i) = {%i};\r\n',y1+1,y1+1);
if Nw==0
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i,%i,%i};\r\n',y1+2,6*(Na-1)+1,6*(Na-1)+2,6*(Na-1)+3,6*(Na-1)+4,6*(Na-1)+5,6*(Na-1)+6);    
%fprintf(fid,'Line Loop(%i) = {%i,%i};\r\n',y1+2,s8+1,s8+2);
fprintf(fid,'Plane Surface(%i) = {%i};\r\n',y1+2,y1+2);
else
 %fprintf(fid,'Line Loop(%i) = {%i,%i};\r\n',y1+2,s8+13,s8+14);   
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i,%i,%i};\r\n',y1+2,s8+3,s8+4,s8+5,s8+8,s8+9,s8+10);
fprintf(fid,'Plane Surface(%i) = {%i};\r\n',y1+2,y1+2);    
end
y2=y1+2;
else
  y2=y1  
end
%Ntaf trailing edge  arası 

%2 kesikli çizgi
for i=1:fT
cc=i*(L0)+1;
cc1=(i-1)*(L0)+1;
    if i==1
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y2+i,s3+i,(cc-1)*6+3,-(s1+i),-(3));
fprintf(fid,'Surface(%i) = {%i};\r\n',y2+i,y2+i);
    elseif i==fT
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y2+i,s3+i,(Na-1)*6+3,-(s1+i),-((cc1-1)*6+3));
fprintf(fid,'Surface(%i) = {%i};\r\n',y2+i,y2+i);     
    else 
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y2+i,s3+i,(cc-1)*6+3,-(s1+i),-((cc1-1)*6+3));
fprintf(fid,'Surface(%i) = {%i};\r\n',y2+i,y2+i);    
    end
end
y3=y2+fT;
for i=1:fT
cc=i*(L0)+1;
cc1=(i-1)*(L0)+1;
    if i==1
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y3+i,s4+i,-((cc-1)*6+4),-(s1+fT+i),(4));
fprintf(fid,'Surface(%i) = {%i};\r\n',y3+i,y3+i);
    elseif i==fT
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y3+i,s4+i,-((Na-1)*6+4),-(s1+fT+i),((cc1-1)*6+4));
fprintf(fid,'Surface(%i) = {%i};\r\n',y3+i,y3+i);                           

    else 
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y3+i,s4+i,-((cc-1)*6+4),-(s1+fT+i),((cc1-1)*6+4));
fprintf(fid,'Surface(%i) = {%i};\r\n',y3+i,y3+i);    
    end
end

y4=y3+fT;
%winglet üst ve alt yüzey

%fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y4+1,-(s7+1),(s12+1),(s7+2),-(s8+1));
%fprintf(fid,'Surface(%i) = {%i};\r\n',y4+1,y4+1);


%fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y4+2,-(s7+1),-(s12+2),(s7+2),(s8+2));
%fprintf(fid,'Surface(%i) = {%i};\r\n',y4+2,y4+2);

%winglet let surface
if Nw==0
y5=y4;
else
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y4+1,s8+1,s8+3,-(s8+6),-(a+6*(Na-1)+1));
fprintf(fid,'Surface(%i) = {%i};\r\n',y4+1,y4+1);
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y4+2,s8+6,s8+4,-(s8+7),-(a+6*(Na-1)+2));
fprintf(fid,'Surface(%i) = {%i};\r\n',y4+2,y4+2);
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y4+3,s8+7,s8+5,-(s8+2),-(a+6*(Na-1)+3));
fprintf(fid,'Surface(%i) = {%i};\r\n',y4+3,y4+3);

fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y4+4,s8+2,s8+8,-(s8+11),-(a+6*(Na-1)+4));
fprintf(fid,'Surface(%i) = {%i};\r\n',y4+4,y4+4);
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y4+5,s8+11,s8+9,-(s8+12),-(a+6*(Na-1)+5));
fprintf(fid,'Surface(%i) = {%i};\r\n',y4+5,y4+5);
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y4+6,s8+12,s8+10,-(s8+1),-(a+6*(Na-1)+6));
fprintf(fid,'Surface(%i) = {%i};\r\n',y4+6,y4+6);
y5=y4+6;
end
%leading edge -ltaf arası up    1
for i=1:fL
cc=i*(L0)+1;
cc1=(i-1)*(L0)+1;
if i==1
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y5+i,s9+i,(cc-1)*6+1,-(s10+i),-(1));
fprintf(fid,'Surface(%i) = {%i};\r\n',y5+i,y5+i);   
elseif i==fL  
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y5+i,s9+i,(Na-1)*6+1,-(s10+i),-((cc1-1)*6+1));
fprintf(fid,'Surface(%i) = {%i};\r\n',y5+i,y5+i);   
else
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y5+i,s9+i,(cc-1)*6+1,-(s10+i),-((cc1-1)*6+1));
fprintf(fid,'Surface(%i) = {%i};\r\n',y5+i,y5+i);
end
end
%%
y6=y5+fL;
%ltaf -Ntaf arası up2
for i=1:fL
cc=i*(L0)+1;
cc1=(i-1)*(L0)+1;
if i==1
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y6+i,s10+i,(cc-1)*6+2,-(s3+i),-(2));
fprintf(fid,'Surface(%i) = {%i};\r\n',y6+i,y6+i);   
elseif i==fL  
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y6+i,s10+i,(Na-1)*6+2,-(s3+i),-((cc1-1)*6+2));
fprintf(fid,'Surface(%i) = {%i};\r\n',y6+i,y6+i);   
else
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y6+i,s10+i,(cc-1)*6+2,-(s3+i),-((cc1-1)*6+2));
fprintf(fid,'Surface(%i) = {%i};\r\n',y6+i,y6+i);
end
end

%%
y7=y6+fL;
%%%%%lead ltaf LOW
%1
for i=1:fL
cc=i*(L0)+1;
cc1=(i-1)*(L0)+1;
if i==1
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y7+i,s9+i,-((cc-1)*6+6),-(s11+i),(6));
fprintf(fid,'Surface(%i) = {%i};\r\n',y7+i,y7+i);   
elseif i==fL  
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y7+i,s9+i,-((Na-1)*6+6),-(s11+i),((cc1-1)*6+6));
fprintf(fid,'Surface(%i) = {%i};\r\n',y7+i,y7+i);   
else
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y7+i,s9+i,-((cc-1)*6+6),-(s11+i),((cc1-1)*6+6));
fprintf(fid,'Surface(%i) = {%i};\r\n',y7+i,y7+i);
end
end
y8=y7+fL;
%ltaf ntaf
for i=1:fL
cc=i*(L0)+1;
cc1=(i-1)*(L0)+1;
if i==1
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y8+i,s11+i,-((cc-1)*6+5),-(s4+i),(5));
fprintf(fid,'Surface(%i) = {%i};\r\n',y8+i,y8+i);   
elseif i==fL  
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y8+i,s11+i,-((Na-1)*6+5),-(s4+i),((cc1-1)*6+5));
fprintf(fid,'Surface(%i) = {%i};\r\n',y8+i,y8+i);   
else
fprintf(fid,'Line Loop(%i) = {%i,%i,%i,%i};\r\n',y8+i,s11+i,-((cc-1)*6+5),-(s4+i),((cc1-1)*6+5));
fprintf(fid,'Surface(%i) = {%i};\r\n',y8+i,y8+i);
end
end
%%

%%
%mesh size
W7=[];
for v=1:fL/2-1
if (v==fL/2-1)
        k1=x11;
    else 
   k1=x1;
end 
 W7=cat(2,W7,k1);
end

rr=(L0*2-1)*4-(Ll*4);

fprintf(fid,'Transfinite Curve {%i:%i}= %i Using Progression 1; \r\n',s1+1,s2,round(200/fL));

fprintf(fid,'Transfinite Curve {%i:%i}= %i Using Progression 1; \r\n',s3+1,s4+fL,round(100/fL));


fprintf(fid,'Transfinite Curve {%i:%i}= %i Using Progression 1; \r\n',s10+1,s11+fL,round(300/fL));


fprintf(fid,'Transfinite Curve {%i:%i}= %i Using Progression 1; \r\n',s9+1,s9+fL,round(2000/fL)); 


fprintf(fid,'Transfinite Curve {%i:%i}= 50 Using Progression 1; \r\n',7,Na*6-6);
fprintf(fid,'Transfinite Curve {%i,%i}= 150 Using Progression 1; \r\n',s7+1,s7+2);

%fprintf(fid,'Transfinite Curve {%i:%i}= 100 Using Progression 1; \r\n',s2+1,s2+round(fL/2));
if Nw==0
ne=1;
else
fprintf(fid,'Transfinite Curve {%i,%i}= %i Using Progression 1; \r\n',s8+13,s8+14,round(200/fL)*6);
fprintf(fid,'Transfinite Curve {%i}= %i Using Progression 1; \r\n',s8+1,Rw*round(2000/fL)*fL); 
fprintf(fid,'Transfinite Curve {%i,%i}= %i Using Progression 1; \r\n',s8+6,s8+12,round(300/fL));
fprintf(fid,'Transfinite Curve {%i,%i}= %i Using Progression 1; \r\n',s8+7,s8+11,round(100/fL));
%fprintf(fid,'Transfinite Curve {%i,%i}= 250 Using Progression 1; \r\n',s12+1,s12+2); 
fprintf(fid,'Transfinite Curve {%i}= %i Using Progression 1; \r\n',s8+2,Rw*round(200/fL)*fL);
end

%fprintf(fid,'Transfinite Curve {%i,%i,%i,%i}= 100 Using Progression 1; \r\n',3,4,(Na-1)*4+4,(Na-1)*4+4);
%fprintf(fid,'Transfinite Curve {%i:%i}= 100 Using Progression 1; \r\n',s5+1,s5+fL); %ntaf ile trailing arası ayıran çizgiler up
%fprintf(fid,'Transfinite Curve {%i:%i}= 100 Using Progression 1; \r\n',s6+1,s6+fL); %ntaf ile trailing arası ayıran çizgiler low


%fprintf(fid,'Transfinite Curve {%i:%i}= 20 Using Progression 1 ;\r\n',s1+1,s1+fL);% slot çizgi çizgi

%fprintf(fid,'Transfinite Curve {%i:%i}= 20 Using Progression 1 ;\r\n',s9+1,s9+fL);% leading edge çizgi

%fprintf(fid,'Transfinite Curve {%i:%i}= 20 Using Progression 1 ;\r\n',s3+1,s3+fL);% ntaf çizgi uıp edge çizgi
%fprintf(fid,'Transfinite Curve {%i:%i}= 20 Using Progression 1 ;\r\n',s4+1,s4+fL);% ntaf çizgi low edge çizgi

%fprintf(fid,'Transfinite Curve {%i,%i}= 10 Using Progression 1 ;\r\n',s7+1,s7+2); %winglet leading and trailing edege
%fprintf(fid,'Transfinite Curve {%i,%i,%i,%i}= 100 Using Progression 1; \r\n',);


%a212=[((L0*2-1)*4)+2:2*rr:((L0*(fL-1)-1)*4)+2];
  
%spaceFormat7=[y12 W7 z12 zx1];
%fprintf(fid,spaceFormat7,a212);

 %a213=[((L0*1-1)*4)+2:2*rr:((L0*(fL-2)-1)*4)+2];
%fprintf(fid,spaceFormat7,a213);

%a214=[((L0*2-1)*4)+5:2*rr:((L0*(fL-1)-1)*4)+5];
%fprintf(fid,spaceFormat7,a214);

 %a215=[((L0*1-1)*4)+5:2*rr:((L0*(fL-2)-1)*4)+5];
%fprintf(fid,spaceFormat7,a215);


%y12='Plane Surface(%i) = {';
%W3=[];
%a3=[1:1:Nf,yyz,yyz+2];
%a4=[Nf:1:2*Nf-2,1,yyz+1,yyz+3];
%for i=1:1:Nf+2
%if (i==Nf+2)
%        k=x11;
%    else 
%   k=x1;
%  
%end 
% W3=cat(2,W3,k);
%end
%spaceFormat=[y12 W3 z1 zx1];
 %fprintf(fid,spaceFormat,yx,a3); 
 %fprintf(fid,spaceFormat,yx+1,a4); 
 
%fprintf(fid,'Surface Loop(%i) = {%i, %i};\r\n',1,1,2);
plot3(ARFLWuy(:,1),ARFLWuy(:,2),ARFLWuy(:,3),'.')
hold on 
plot3(ARFLWly(:,1),ARFLWly(:,2),ARFLWly(:,3),'.')
