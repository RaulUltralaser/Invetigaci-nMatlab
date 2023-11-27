function [us] = MeasureData(t)
%------------------------------------------------------------------
% Purpose:
% Returns the state of the unknown DPS for all nodes
%
% Sinopsys:
% us = MeasureData(t)
%
% Variable description:
% t - current simulation time
%
%------------------------------------------------------------------
persistent x y z
if isempty(x)
	p = 10;
	X = linspace(0,1,p);
	Y = linspace(0,1,p);
	Z = linspace(0,1,p);

	[x,y,z] = meshgrid(X,Y,Z);
end
us(:,:,:) = 0.5*sin(pi*x).*sin(pi*y).*sin( sin(pi*sqrt(2)*t)/(pi*sqrt(2)) ... 
    	                              -cos(2*pi*z).*sin(pi*sqrt(6)*t)/(pi*sqrt(6)));
us = reshape(us,[],1);