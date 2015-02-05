classdef LTIADPsys < handle
	
	properties
		A
		B
		C
		R
		K0 % initial gain
		K  % current gain
	end
	
	methods
		function obj = LTIADPsys(A,B,C,R,K0)
			% To construct an LTIADPsys object, use the syntax:
			% obj = LTIADPsys(A,B,C,R,K0)
			% where
			% (A,B) is a pair of stabilizable matrices
			% C, with (A,C) observable is to create Q such that Q = C'*C.
			% R is the weight for the control input
			% K0 is the initial stabilizing feedback gain
			
			% Validate the input paramters
			if size(A,1) ~= size(B,1)
				error('Matrices A and B must have the same number of rows');
			end
			
			obj.A = A;
			obj.B = B;
			obj.C = C;
			obj.R = R;
			obj.K0 = K0;
			obj.K = K0;
		end
		
		function dx = odefcn(obj,t,x,flag)
			% flag:  1. exploration (learning phase), noise on
			%		 0. exploitation (implementation phase), noise off
			
			u = obj.K*x + explorationNoise(t, flag);
			dx = obj.A*x - obj.B*u;
		end
		
		function [t,y] = simulate(obj,x0,t0,tf,flag)
			% This method simulats the given system on a specified interval
			% [t0,tf], with the initial condition x0. The argument flag is
			% to control if the noise is applied or not.
			[t,y] = ode45(@(t,x) obj.odefcn(t,x,flag), [0,10], 1);
			if nargout < 2
				figure
				plot(t,y)
			end
		end
		
		function update(obj)
			% This method updates once the control gain
			% blablabla
		end
	end
	
	
end


%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = explorationNoise(t, flag)
% flag:  1. exploration (learning phase), noise on
%		 0. exploitation (implementation phase), noise off

if flag == 0
	e = 0;
	return
else
	% use the noise below, or custmize your own noise generator
	e  = sin(t);
end

end
