%finite element analysis for the euler bernouli beam

clear A Kc Mc Dc nq Mbc Kbc y yo v
close all
ne=21;           %numero de elementos                  
nne=2;          %numero de nodos por elemento
dof=2;          %Grado de libertad -> Nùmero de variables por nodo
L=1;            %Longitud de la viga
Le=L/ne;        %H distancia entre elementos
E=69E7;        %Nm %Modulo de Young          -> Caracterizacón mecánica
I=1.66667E-11;            %Momento de inercia
area=5E-5;         %area sección transversal  m^2
rho=2710;          %densidad del metal        kg/m^3
nn=ne+1;        %Número de nodos
nq=dof*nn;      %Número de variables totales

%Obtenemos alfa y beta con caracterización mecánica 
alpha=0.01;
beta=0.0001;

x=0:Le:L;
x=x';
yo=zeros(nn,1);
nodes=zeros(ne,nne);
for i=1:ne
nodes(i,:)=[i,i+1];
end

KK=zeros(nq,nq);
MM=zeros(nq,nq);

imass=1;                            %consistent mass matrix

for e=1:ne
[K,Mg]=febeam1(E,I,Le,area,rho,imass);
index=feeldof(nodes(e,:),nne,dof);
KK=feasmbl1(KK,K,index);            %MATRIZ DE TODO EL SISTEMA
MM=feasmbl1(MM,Mg,index);
end


%%boundary conditions
%simply supported beam
bc=1;
nbc=length(bc);
bcval=zeros(1,dof*nbc);
ibc=feeldof(bc,nbc,dof);


%inputs and outputs
f1=nn;
% f2=nn-1;
h=1:nn;
% h2=4;
%traction vector
qi=f1;  %traction surface NODES
qo=h;  %measurable nodes
m=length(qi);
p=length(qo);

fv=-0.1;  %displacement force 
fa=0;  %angle force


yr=-0.1;
%yr=[-0.1;-0.3];



%phi=[fv,fa,fv,fa]';
phi=[-0.1;0];


iqi=feeldof(qi,m,dof);
iqo=feeldof(qo,p,dof);

F=zeros(nq,1);
F(iqi)=phi;

FF=zeros(nq,m);
fi=[1,0]';
i=1;
for j=1:m
    FF(iqi(i:i+1),j)=fi;
    i=i+dof;
end

H=zeros(p,nq);
hi=[1,0];
i=1;
for j=1:p
H(j,iqo(i:i+1))=hi;
i=i+dof;
end



%static analysis
[Ks,F]=feaplyc2(KK,F,ibc,bcval);
q=Ks\F;         % X = A/B = A^-1*B     || X = A/B = A*B^-1

%arrangement of displacement vector
j=1;
for i=1:nn
v(i,1)=q(j);
angle(i,1)=q(j+1);
j=j+2;
end
G=1;


%display results
y=yo+G*v;

figure()
plot(x,y,'k');
hold on
plot(x(f1),yo(f1),'or')  %force nodes in red
plot(x(h),yo(h),'ob')  %output nodes in blue
%plot(XF([100,103:104]),YF([100,103:104]),'*r')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%dynamic analysis
[KK,MM]=feaplycs(KK,MM,ibc);
nq=dof*(nn-nbc);
Hbc=zeros(p,nq);
Fbc=zeros(nq,m);
Kbc=zeros(nq,nq);
Mbc=zeros(nq,nq);


%miss remove bc!!
%REMOVE fixed nodes (Dirichlet boundary conditions)
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

Dbc=alpha*Mbc+beta*Kbc;


% ESPACIO DE ESTADOS
%system with bc 
A=[zeros(nq,nq),eye(nq);
   -Mbc\Kbc  ,-Mbc\Dbc];

B=[zeros(nq,m);
    Mbc\Fbc];
C=[Hbc,zeros(p,nq)];
Cq=diag([ones(1,nq),zeros(1,nq)]);


l=eig(A);


%%
% for i = 1:length(Cq)
%     for j=1:length(Cq)
%         if mod(Cq(i,j),2) ~= 0 && i<=8
%             Cq_new(i,j) = 1;
%         elseif mod(i,2)==0 && i<=8
%             Cq_new(i,j) = 0;
%         else 
%             Cq_new(i,j) = 0;
%         end
%     end
% end
% 
% 
% for i = 1:length(Cq)
%     for j=1:length(Cq)
%         if i==j && mod(i,2)==0
%              Cq_new(i,j) = 0;
%         end
%     end
% end