function [du, dW, dV] = DNN(u, us, W1, V1, h)
	global sigmoid K1 K2 P V0 l Lambda A
	Delta   = u - us;
	sigma   = sigmoid(1,V1*u);
	D_sigma = diag(sigma.*(1-sigma));
	%
	dW      = h*(-K1*P*Delta*sigma' + K1*P*Delta*u'*(V1-V0)'*D_sigma);
	dV		= h*(-K2*D_sigma'*W1'*P*Delta*u' - l/2*K2*Lambda*(V1-V0)*(u*u'));
	du		= h*(A*u + W1*sigma);
end