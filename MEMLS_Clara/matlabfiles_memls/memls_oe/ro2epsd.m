   function [epsi,epsii] = ro2epsd(roi,Ti,freq)

%   calculates the dielectric permittivity from 
%   density for dry snow.
%
%   [epsi,epsii] = ro2epsd(roi,Ti,freq)
%       epsi:  real part of dielectric permittivity
%       epsii: imaginary part of dielectric permittivity
%       roi:   density
%       Ti:    snow temperature in Kelvin
%       freq:  frequency
%
%   Version history:
%      1.0    wi 15.7.95
%      2.0    wi 12.11.97  enhanced with Polder and van Santen Equations (see Polder.m)
%
%   Uses:
%       epsice, epsr, polder
%
%
%
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

eice = epsice(Ti,freq);

epsi = epsr(roi);

% imaginary part after Tiuri 84
%epsii = eice.*(0.52.*roi + 0.62.*(roi.^2));

% imaginary part after Polder and van Santen 1946 (Effective-Medium Approx)
f = roi ./ 0.917;
ei = 3.185;
N = max(size(roi));
A = zeros(N,1);
epsp = zeros(N,1);
A = A + 0.3;
for i=1:N
   if f(i) < 0.55
      A(i) = 0.476 - 0.64 * f(i);
   end
   if f(i) <= 0.333
      A(i) = 0.1 + 0.5 * f(i);
   end
   %epsp(i) = polder(A(i),ei,f(i));
end
epsp = epsi;
A3 = 1 - 2 .* A;
ea = (epsp .* (1-A)) + A;
ea3 = epsp .* (1-A3) + A3;
K1 = (ea  ./ (ea+A   .* (ei-1))).^2;
K3 = (ea3 ./ (ea3+A3 .* (ei-1))).^2;
Ksq = (2 .* K1 + K3) ./ 3;
epsii = sqrt(epsi) .* eice .* Ksq .* f; 
