%---------------------------------------------------------------------------
% FEM_DNN_tetrahedral_element.m
% Main program for DPS identification using FEM and DNN
%
%
%---------------------------------------------------------------------------
clear all
close all
clc
addpath('.\ejemplos')
%--------------------------------------
% input data for control parameters
%--------------------------------------
npa   = 10;						% nodes per axis
nel   = 6*(npa-1)^3;    		% number of elements
nnel  = 4;						% nodes per element
ndof  = 1; 						% dofs per node
nnode = npa^3;					% number of nodes in system
sdofs = nnode*ndof;				% total system dofs
%
%----------------------------------------------
% Create the nodes and give them their coordinates
%----------------------------------------------
% generateMesh()
%CreateMesh
%
%----------------------------------------
% Simulation parameters
%----------------------------------------
ts = 2;									% simulation time
h  = 0.01;								% sample time
tk = 0:h:ts;							% time vector
N  = length(tk);						% number of iterations
%
%------------------------------------
% initialization of matrices and vectors
%------------------------------------
ff    = zeros(sdofs,1);				% system vector
fn    = zeros(sdofs,1);				% effective system vector
fsol  = zeros(sdofs,1);             % solution vector
sol   = zeros(1,N+1);               % time-history solution
kk    = zeros(sdofs,sdofs);			% system matrix
mm    = zeros(sdofs,sdofs);			% initialization of system matrix
kn    = zeros(sdofs,sdofs);			% effective system matrix
index = zeros(nnel*ndof,1);         % index vector
%
%------------------------------------------------
% computation of element matrices and vectors and their assembly
%------------------------------------------------
for ielement = 1:nel 				% loop for the total number of elements
	%
	nd = nodes(ielement,:);			% connected nodes for (iel)-th element
	x  = nodal_coord(nd,1)';
	y  = nodal_coord(nd,2)';
	z  = nodal_coord(nd,3)';		% coordinates of the 4 nodes
	%
	index = feeldof(nd,nnel,ndof);		% extract system dofs for the element
	%
	k = felp3dt4(x,y,z);				% compute element matrix
    m = felpt3t4(x,y,z);                % compute element matrix
	%
	kk = feasmbl1(kk,k,index);			% assemble element marices
    mm = feasmbl1(mm,m,index);          % assemble element matrices
	%
end
clear ielement k m nd x y z
%%
clc
close all
%----------------------------------------------
% Parameters initialization for DNN
%----------------------------------------------
global V1 W1 sigmoid K1 K2 P V0 l Lambda A 
V1 = 2*rand(nnode,nnode)-1;			% weigth matrix
W1 = 2*rand(nnode,nnode)-1;			% weigth matrix
us = MeasureData(0);                % "real" measurement
u  = us;                            % first state of the system
%
Kmask 			= kk;
% Kmask(Kmask~=0) = 1;
V1		        = V1.*Kmask;
W1              = W1.*Kmask;
V0              = V1;
%
sigmoid = @(b,x)( 1./(1+exp(-b*x)));
I       = eye(nnode);
K1		= 2.2802;
K2		= 2.7468;
l 		= 1.1620;
P       = I;%SPDmatrix(nnode);
Lambda	= SPDmatrix(nnode);
aa      = -25;%-51.1440;
A       = aa*eye(nnode);
%
Q0  = I;%SPDmatrix(nnode);
Q   = Q0 +Lambda;
bW1 = W1*inv(Lambda)*W1';
R   = 2*bW1 +inv(Lambda);

a = [-P*A-A'*P-Q P; P inv(R)];
try chol(a);
    disp("All parameters P, A, Q and R satisfy the conditions ...")
catch ME
    disp('Matrix is not symmetric positive definite')
end
%[-P*A-A'*P-Q P; P inv(R)]
%A = eye(nnode);
%% 
%----------------------------------------
% Algorithm 1: DPS Identification 
%----------------------------------------
measured_u = zeros(npa^3,N);        % real state
approx_u   = zeros(npa^3,N);        % approx state
%
X = linspace(0,1,npa);
Z = linspace(0,1,npa);
[x,z]  = meshgrid(X,Z);
graph1 = zeros(npa,npa);
graph2 = zeros(npa,npa);
error  = zeros(1,N);
figure(1)
tic
for i = 1:N
	W1 = W1.*Kmask;
	V1 = V1.*Kmask;
    %
    %DoPSO
	%
	[k1u, k1W, k1V] = DNN(u,us,W1,V1,h);
	%
	us = MeasureData(tk(i)+0.5*h);
	%
	[k2u, k2W, k2V] = DNN(u+.5*k1u,us,W1+.5*k1W,V1+.5*k1V,h);
	[k3u, k3W, k3V] = DNN(u+.5*k2u,us,W1+.5*k2W,V1+.5*k2V,h);
	%
	us = MeasureData(tk(i)+h);
	%
	[k4u, k4W, k4V] = DNN(u+k3u,us,W1+k3W,V1+k3V,h);
	%
	u  =  u + 1/6*k1u + 1/3*k2u + 1/3*k3u + 1/6*k4u;
	V1 = V1 + 1/6*k1V + 1/3*k2V + 1/3*k3V + 1/6*k4V;
	W1 = W1 + 1/6*k1W + 1/3*k2W + 1/3*k3W + 1/6*k4W;
    
    % Real model
	measured_u(:,i) = us;
	% Neural Network model
    approx_u(:,i)   = u;
    %
    new_u = reshape(measured_u(:,i),npa,npa,npa);
    aprox = reshape(approx_u(:,i),npa,npa,npa);
    %
    for kk = 1:length(Z)
        graph1(:,kk) = diag(new_u(:,:,kk));			% j -> 131
        graph2(:,kk) = diag(aprox(:,:,kk));			% j -> 131
    end
    subplot(221)
    contourf(x,z,graph1,10);
    colormap(jet);
    pcolor(x,z,graph1);
    cb1 = colorbar;
    shading flat
    title('Measured data')
    subplot(222)
    contourf(x,z,graph2,10);
    colormap(jet);
    pcolor(x,z,graph2);
    cb2 = colorbar;
    shading flat
    title('Approximated data')
    %
    subplot(2,2,3)
    error(i) = 1/i*(new_u(5,5,5)-aprox(5,5,5))^2;
    plot(tk(1:i), error(1:i), 'r')
    title(strcat('MSE(t=',num2str(tk(i)),') = ',num2str(error(i))))
    xlabel('time [sec]')
    set(gcf, 'Position',  [300, 100, 650, 500])
    grid on
    %
    subplot(2,2,4)
    error(i) = 1/i*sum((new_u-aprox).^2,'all');
    plot(tk(1:i), error(1:i), 'r')
    title(strcat('MSE(t=',num2str(tk(i)),') = ',num2str(error(i))))
    xlabel('time [sec]')
    set(gcf, 'Position',  [300, 100, 650, 500])
    grid on
    
    pause(0.05)
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
        imwrite(imind,cm,'testnew51.gif','gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,'testnew51.gif','gif','WriteMode','append');
    end
end
toc
% Distributed Parameter System Identification using 
% Finite-Element Differential Neural Networks 
% (O. A. Aguilar, R. Q. Fuentes) pg. 9
%-------------------------------------------------------------------
%--------------------------------------------
% System simulation Parameters
%--------------------------------------------
% ts = 2;                     % simulation time
% t0 = 0;                     % initial time
% T  = 0.01;                  % sample period
% t  = 0:T:ts;                % time vector
% N  = length(t);             % number of iterations
%

%approx_u   = zeros(npa^3,N);        % approx state
%
%W1 = W1.*Kmask;
%V1 = V1.*Kmask;
%
% for i = 1:N-1
% 	% Real model
% 	measured_u(:,i+1) = MeasureData(t(i));
% 	% Neural Network model
%     approx_u(:,i+1)   = sys(@DNN_model, approx_u(:,i), t0+(i-1)*T, T);
% end
% redimensioned_u = reshape(responce,npa,npa,npa);
%% Plots
% X = linspace(0,1,npa);
% Y = linspace(0,1,npa);
% Z = linspace(0,1,npa);
% graph1 = zeros(npa,npa);
% graph2 = zeros(npa,npa);
% [x,z] = meshgrid(X,Z);
% for j = 1:length(t)
%     new_u = reshape(measured_u(:,j),npa,npa,npa);
%     aprox = reshape(approx_u(:,j),npa,npa,npa);
%     for i = 1:length(Z)
%         graph1(:,i) = diag(new_u(:,:,i));			% j -> 131
%         graph2(:,i) = diag(aprox(:,:,i));			% j -> 131
%     end
%     subplot(121)
%     contourf(x,z,graph1,10);
%     colormap(jet);
%     pcolor(x,z,graph1);
%     cb1 = colorbar;
%     shading flat
%     subplot(122)
%     contourf(x,z,graph2,10);
%     colormap(jet);
%     pcolor(x,z,graph2);
%     cb2 = colorbar;
%     shading flat
%     pause(0.5)
% end
%
function x = sys(model, X, t0, T)
    opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
    X    = X';
    %[t, x_] = ode45(model, linspace(t0, t0+T, 3), X, opts);
    [t, x_] = ode45(model, [t0, t0+0.5*T, t0+T], X, opts);
    x       = x_(end,:)';
end
%
function ut = DNN_model(t,u)
    global A sigmoid V1 W1
    ut = A*u + W1*sigmoid(1,V1*u);
end