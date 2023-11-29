% clear all
close all
clear all
clc
%% Parameters
Np     = 500;                             % simulation time (S)
X      = zeros(2, Np);                    % states over time 
hX     = X;                               % approximate state
X(:,1) = [10; -10];                       % initial state
t      = 1:Np;                            % simulation time vector
fs     = 1/100;                           % frequency for the input signals
ut     = [sin(2*pi*t*fs); ...
         -sawtooth(2*pi*t*fs)];           % input signals
t0     = 0;                               % initial time
Ts     = 1;                               % sampling time

% DNN parameters
global W1t W2t sigma phi
hX(:,1) = [-5; -5];
I       = eye(2);
K2      = 10*I;
K1      = 10*I;
P       = [0.06 0.04; 0.04 0.106];
Q0      = I;
R       = [8 2; 2 8];
Q       = [2 1; 1 2];

% Sigmoid functions
sigmoid = @(b,x)( 1./(1+exp(-b*x)));
sigma   = @(xi)( 2*sigmoid(2,xi) - 0.5);
phi     = @(xi)(diag(0.2./(1+exp(-0.2*xi))-0.05));

% Dead-zone function
bmu        = 1;
beta_phi   = 2;
beta_sigma = 2;
mu         = (beta_sigma + beta_phi*bmu)/min(eigs((P^(-1/2)*Q0*P^(-1/2))));
q          = 1;


%% Single layer network
W1t = [1 10; 10 1];
W2t = [0.1 0; 0 0.1];

for k = 1 : Np-1
    % Real model
    X(:,k+1)  = sys(@Real_model, X(:,k), ut(:,k), t0+(k-1)*Ts, Ts); 

    % Neural Network model
    hX(:,k+1) = sys(@NNSL_model, hX(:,k), ut(:,k), t0+(k-1)*Ts, Ts);

    % Error
    delta     = hX(:,k+1) - X(:,k+1);

    % Dead-zone function
%     st = 1 - mu/norm(P^(1/2)*delta);
%     if st < 0; st = 0; end
    st = 1;

    % Learning law
    dW1t = -st*K1*P*delta*sigma(hX(:,k+1))';
    dW2t = -st*K2*P*delta*ut(:,k)'*phi(hX(:,k+1));

    % Weigths updating
    W1t = W1t + dW1t;
    W2t = W2t + dW2t;
end

subplot(1,2,1)
plot([2:k], X(1,2:k), 'b', [2:k], X(2,2:k), 'r'); hold on   % Real model
plot([2:k], hX(1,2:k), 'b--', [2:k], hX(2,2:k), 'r--')      % NN model
set(legend,'Interpreter','latex')
legend('$x_1$', '$x_2$', '$\hat{x}_1$', '$\hat{x}_2$', 'FontSize', 14)
title("Single layer network", 'FontSize', 14)
axis([0 500 0.1 0.8])
xlabel("time [s]", 'FontSize', 14)
ylabel("States", 'FontSize', 14)

global V1t V2t
%% Multilayer network
W1t = [1 1 2; 1 2 1];
W2t = [1 1 2; 1 2 1];
V10 = [1 1 2; 1 2 1]';
V20 = [1 1 2; 1 2 1]';
V1t = W1t';
V2t = V1t;
K3  = eye(3);
K4  = K3;

for k = 1 : Np-1
    X(:,k+1)  = sys(@Real_model,  X(:,k), ut(:,k), t0+(k-1)*Ts, Ts);
    hX(:,k+1) = sys(@NNML_model, hX(:,k), ut(:,k), t0+(k-1)*Ts, Ts);
    delta     = hX(:,k+1) - X(:,k+1);           
    st        = 1;

    % Updating law
    sig     = sigma(V1t*hX(:,k+1));
    D_sigma = diag(2*2*sig.*(1-sig));

    fi      = phi(V2t*hX(:,k+1));
    D_phi   = 0.2*0.2*fi.*(1-fi);

    suma    = [ut(:,k); 0]*D_phi(:,1)' + ...
              [ut(:,k); 0]*D_phi(:,2)' + ...
              [ut(:,k); 0]*D_phi(:,3)';

    l1      = 2;                % l1, l2 > 0
    l2      = 3;
    Lambda1 = eye(3);           % Lambda are positive defined matrices
    Lambda2 = eye(3);
    
    dW1t = -st*K1*P*delta*sig' ...
           +st*K1*P*delta*hX(:,k+1)'*(V1t-V10)'*D_sigma;
    dW2t = -st*K2*P*delta*(fi*[ut(:,k); 0])' ...
           +st*K2*P*delta*hX(:,k+1)'*(V2t-V20)'*suma;
    dV1t = -K3*D_sigma'*W1t'*P*delta*hX(:,k+1)' ...
           -(l1/2)*K3*Lambda1*(V1t-V10)*hX(:,k+1)*hX(:,k+1)';
    dV2t = -K4*suma*W2t'*P*delta*hX(:,k+1)' ...
           -(q*l2*bmu/2)*K4*Lambda2*(V2t-V20)*hX(:,k+1)*hX(:,k+1)';

    W1t = W1t + dW1t;
    W2t = W2t + dW2t;
    V1t = V1t + dV1t;
    V2t = V2t + dV2t;
end
subplot(1,2,2)
plot([2:k], X(1,2:k), 'b', [2:k], X(2,2:k), 'r'); hold on   % Real model
plot([2:k], hX(1,2:k), 'b--', [2:k], hX(2,2:k), 'r--')      % NN model
set(legend,'Interpreter','latex')
legend('$x_1$', '$x_2$', '$\hat{x}_1$', '$\hat{x}_2$', 'FontSize', 14)
title("Multilayer network", 'FontSize', 14)
xlabel("time [s]", 'FontSize', 14)
ylabel("States", 'FontSize', 14)
axis([0 500 0.1 0.8])
%sgtitle("Example 2.1 (Poznyak)", 'FontSize', 16)
set(gcf, 'Position',  [300, 100, 1050, 400])
%%
function x = sys(model, X, u, t0, T)
    opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
    X = X';
    [t, x_] = ode45(model, linspace(t0, t0+T, 20), X, opts, u);
    x       = x_(end,:)';
end
function X = Real_model(t, x, u)
    x1 = x(1);
    x2 = x(2);
    u1 = u(1);
    u2 = u(2);
    
    dx1 = -5*x1 +3*sign(x2) +u1;
    dx2 = -10*x2 +2*sign(x1) +u2;
    X   = [dx1; dx2]; 
end
function dhX = NNSL_model(t,hxt,ut)
    global W1t W2t sigma phi
    A   = [-15 0; 0 -10];
    dhX = A*hxt + W1t*sigma(hxt) + W2t*phi(hxt)*ut(:,1);
end
function dhx = NNML_model(t,hxt,ut)
    global V1t V2t W1t W2t sigma phi
    A   = [-15 0; 0 -15];
    dhx = A*hxt + W1t*sigma(V1t*hxt) + W2t*phi(V2t*hxt)*[ut; 0];
end