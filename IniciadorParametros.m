%clear Kc Mc Dc mu mu2 betap betav
clearvars
close all


f=20;        %que nodo quiero leer
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
I=a^4/12;   %m4 moment of inertia
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


%% DNN
%----------------------------------------------
% Parameters initialization for DNN
% %----------------------------------------------
% global V1 W1 sigmoid K1 K2 P V0 l Lambda A 
% Agregue NN A TODAS LAS VARIABLES PARA EVITAR CONFUSIONES DE ESTA PARTE Y
% LA ANTERIOR
nnode=80;                              %Está relacionado a FF
V1NN = 2*rand(nnode,nnode)-1;			      % weigth matrix
W1NN = 2*rand(nnode,nnode)-1;			      % weigth matrix
% us = MeasureData(0);                % "real" measurement
% u  = us;                            % first state of the system
%
KmaskNN 			= [zeros(nq,nq),eye(nq);
                        Kbc  ,Dbc]; ;
% Kmask(Kmask~=0) = 1;
V1NN		        = V1NN.*KmaskNN;
W1NN                = W1NN.*KmaskNN;
V0NN                = V1NN;
W0NN                = W1NN;                     %%YO AGREGUE ESTO PARA INICIALIZAR EN SIMULINK
%
% sigmoid = @(b,x)( 1./(1+exp(-b*x)));
hNN       =0.01;%sample time   cambie el nombre originalmente (h)
K1NN	  = 2.2802;
K2NN	  = 2.7468;
INN       = eye(nnode);
lNN 	  = 1.1620;
PNN       = INN;%SPDmatrix(nnode);
LambdaNN  = SPDmatrix(nnode);
aaNN      = -25;%-51.1440;
ANN       = aaNN*eye(nnode);   


function A = SPDmatrix(size)
    % Generate a random symmetric matrix
    A = randn(size, size);

    % Make the matrix symmetric
    A = 0.5 * (A + A');

    % Make the matrix positive definite
    A = A + size * eye(size);
end


