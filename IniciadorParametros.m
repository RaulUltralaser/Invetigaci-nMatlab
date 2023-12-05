clear Kc Mc Dc mu mu2 betap betav
close all

% Valores de la viga
ne=21;      %Nodos
nne=2;      %Nodos Elementos
dof=2;      %Grados de libertad
L=1;        %Large %m
Le=L/ne;    %Element distance
E=69e6;     %Young's module N/m2

a=0.0015;   
d=2e-3;     %tickness m
b=25e-3;    %broadness m

area=b*d;   %cross-sectional area
%area=a^2;
I=a^4/12;   %m4 moment of area
%Iy=b*d^3/12; %area moment for forces on z-axis Fz

%Iz=b^3*d/12;
%I=Iy;
rho=2710;    %density Kg/m3
imass=1;     %consistent mass matrix
nn=ne+1;
N=nn*dof;
Rm=1e-3;
Rk=1e-4;

x=0:Le:L;
x=x';
yo=zeros(nn,1);
nodes=zeros(ne,nne);
for i=1:ne
nodes(i,:)=[i,i+1];
end

%Initialize stiffness and mass matrix
KK=zeros(N,N);
MM=zeros(N,N);

%Arrange the matrices
for e=1:ne
    [K,M]=febeam1(E,I,Le,area,rho,imass);
    index=feeldof(nodes(e,:),nne,dof);
    KK=feasmbl1(KK,K,index);
    MM=feasmbl1(MM,M,index);
end



%simply supported beam
bc=[1,nn]; 
nbc=length(bc);
bcval=zeros(1,dof*nbc);
ibc=feeldof(bc,nbc,dof);
%traction vector

qi=2:nn-1; %controled nodes 
qo=2:nn-1;  %measurable nodes  

m=length(qi);
p=length(qo);

%fd=1;  %displacement force 
%fa=0;  %angle force


iqi=feeldof(qi,m,dof);
iqo=feeldof(qo,p,dof);

FF=zeros(N,m);

fi=[1,0]';  %displacement force
%fi=[0,1]'; %angular force
i=1;
for j=1:m
    FF(iqi(i:i+1),j)=fi;
    i=i+dof;
end

H=zeros(p,N);
hi=[1,0];
i=1;
for j=1:p
H(j,iqo(i:i+1))=hi;
i=i+dof;
end

%%dyanmic analysis
[KK,MM]=feaplycs(KK,MM,ibc);

%%%%%%%%%%%%%%
%%Remove fixed nodes (Dirichlet boundary conditions)
nq=dof*(nn-nbc);
Hbc=zeros(p,nq);
Fbc=zeros(nq,m);

Kbc=zeros(nq,nq);
Mbc=zeros(nq,nq);
qq=1:2*nn;
qq(ibc)=[];


for i=1:nq
    for j=1:nq
        Kbc(i,j)=KK(qq(i),qq(j));
        Mbc(i,j)=MM(qq(i),qq(j));
    end
    Fbc(i,:)=FF(qq(i),:);
    Hbc(:,i)=H(:,qq(i));
end

Dbc=Rm*Mbc+Rk*Kbc;


%system with bc 
A=[zeros(nq,nq),eye(nq);
   -Mbc\Kbc  ,-Mbc\Dbc];   
B=[zeros(nq,m);
    Mbc\Fbc];


%Fbcw2=[1,0,1,0,1,0,1,0,1,0,1,0,1,0]';
%Bw2=[zeros(nq,1);Mbc\Fbcw2];  %Usar este vector para H_infty

C=[Hbc,zeros(p,nq)];
C2=[zeros(p,nq),Hbc];


[Phi,W2]=eig(Kbc,Mbc);
Z=Rm*eye(nq)+Rk*W2;
z=diag(Z);
w2=diag(W2);
w=real(sqrt(w2));
zeta=z./(2*w);
Psi2=[Phi'*Mbc,zeros(nq);zeros(nq),Phi'*Mbc];
Phi2=[Phi,zeros(nq);zeros(nq),Phi];
iPhi=inv(Phi);

% %modal reduction
r=5;
k=2*r;
W2r=W2(1:r,1:r);
Zr=Z(1:r,1:r);
Phir=Phi(:,1:r);
iPhir=iPhi(1:r,:);
Phir2=[Phir,zeros(nq,r);zeros(nq,r),Phir];
iPhir2=[iPhir,zeros(r,nq);zeros(r,nq),iPhir];

Am2=[zeros(r,r),eye(r);
    -W2r, -Zr];
Bm2=[zeros(r,m);                    %%%This will be useful for optimal SAP
    Phir'*Fbc];

Cm2=[Hbc*Phir,zeros(p,r)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[zeta_min,i]=min(zeta);
cr=zeros(1,r);
i=1;  %mode
cr(i)=1;
Cv=[zeros(1,r),cr];
Cp=[cr,zeros(1,r)];

Cm=[ones(1,r),zeros(1,r)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Sensor Actuator Performance Measure (SAP)

zetad=1/sqrt(2);
zetad=0.7;
zetad=1;
zetag=zetad-zeta(i);



PhiB=Phir'*Fbc;
PhiC=Hbc*Phir;
PB=Phi'*Fbc;
PC=Hbc*Phi;

L1=max(abs(PhiB(i,:)));
L2=min(abs(PhiB(i,:)));

rmodes=1:r;
rmodes(i)=[]; %residual modes


mu0=zeros(p,m);
mu1=zeros(p,m);
beta=zeros(p,m);
betav=zeros(p,m);
au=zeros(1,r-1);


for h=1:p
    ci=abs(PhiC(h,i));
    
    for f=1:m
        bi=abs(PhiB(i,f));
        if (ci*bi)<1e-9 
            mu0(h,f)=inf;     
        else
            for j=1:r-1
                rm=rmodes(j);
                au(j)=abs(PhiC(h,rm))*abs(PhiB(rm,f))/abs(w2(rm)-w2(i));
                c2(j) =abs(PhiC(h,rm))*abs(PhiB(rm,f))/sqrt((w2(rm)-w2(i))^2+(2*zeta(rm)*w(rm)*w(i))^2);
                c8(j)=abs(PhiC(h,rm))*abs(PhiB(rm,f))/(2*zeta(rm)*w(rm));           
            end      
        end
        mu0(h,f)=L1*ci/(2*zeta(i)*w2(i));
        mu1(h,f)=L1*sum(au)/bi;   
        beta(h,f)=L1*ci/(2*zetad*w2(i))+mu1(h,f);
        betav(h,f)=L1*ci/(2*zetad*w(i))+mu1(h,f);
    end
end

%mu  %Utilizamos |G^j|, que es menor que la norma infinita de G
%mu2 %Utlizamos norma infinita

%% DNN
%----------------------------------------------
% Parameters initialization for DNN
%----------------------------------------------
% global V1 W1 sigmoid K1 K2 P V0 l Lambda A 
nnode=44;                              %NO SÉ QUE SE SUPONE QUE SEA NNODE AQUÍ
V1 = 2*rand(nnode,nnode)-1;			      % weigth matrix
W1 = 2*rand(nnode,nnode)-1;			      % weigth matrix
% us = MeasureData(0);                % "real" measurement
% u  = us;                            % first state of the system
%
Kmask 			= KK;
% Kmask(Kmask~=0) = 1;
V1		        = V1.*Kmask;
W1              = W1.*Kmask;
V0              = V1;
W0              =W1;                     %%YO AGREGUE ESTO PARA INICIALIZAR EN SIMULINK
%
% sigmoid = @(b,x)( 1./(1+exp(-b*x)));
h=0.01;                       %sample time   
K1		= 2.2802;
K2		= 2.7468;
I       = eye(nnode);
l 		= 1.1620;
P       = I;%SPDmatrix(nnode);
Lambda	= SPDmatrix(nnode);
l 		= 1.1620;


function A = SPDmatrix(size)
    % Generate a random symmetric matrix
    A = randn(size, size);

    % Make the matrix symmetric
    A = 0.5 * (A + A');

    % Make the matrix positive definite
    A = A + size * eye(size);
end


