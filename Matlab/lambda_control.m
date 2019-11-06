function [validity]= lambda_control(soil_parameters)
%Check if K(h) is strictly >0
%   This program controls the validity of the soil parameter set for each node.
%   The condition is that the derivative of the hydraulic conductivity
%   as a function of relative water content must be all over positive.
%   This condition will be all over satisfied if and only if it is satisfied for high
%   pressure heads.
%IN:
%   soil_parameters
%OUT:
%  validity
%CALLS:
%   none
%CALLED BY:
%   wavemat.m
%----------------------------
%by Lambot S.,2000

validity = 1;
h = 10000;
for node = 1:size(soil_parameters,1)
    wcr    = soil_parameters(node,1);
    wcs    = soil_parameters(node,2);
    alpha  = soil_parameters(node,3);
    n      = soil_parameters(node,4);
    ks     = soil_parameters(node,5);
    lambda = soil_parameters(node,6);
   
    m  = 1-1/n;
    Se = (1+(alpha*h)^n)^(-m);
    a  = (1-(1-Se^(1/m))^m)^2;
    b  = 2*(1-(1-Se^(1/m))^m)*(-m*(1-Se^(1/m))^(m-1))*(-1/m*Se^(1/m-1));
  
   if (a~=0)&(b~=0)
      if (-b*Se/a)>lambda
         validity = 0;
      end
   end
end
