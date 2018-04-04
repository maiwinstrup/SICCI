function [gbih,gbiv,gs6,ga2i] = sccoeff(roi,Ti,pci,freq,Wi,gai,sccho)

%   calculates the scattering coefficient from structural parameters
%     different algorithms can be chosen, by changing "sccho"
%
%   [gbih,gbiv,gs6,ga2i] = sccoeff(roi,Ti,pci,freq,Wi,gai,sccho)
%       gbih:  2-flux scattering coefficient at h pol
%       gbiv:  2-flux scattering coefficient at v pol
%       gs6:   6-flux scattering coefficient
%       ga2i:  2-flux absorption coefficient
%       roi:   density
%       Ti:    physical temperature
%       pci:   correlation length
%       freq:  frequency
%       Wi:    wetness
%       gai:   absorption coefficient
%       sccho: scattering coefficient algorithm chosen
%
%   Version history:
%      1.0b    wi 15.7.95
%      1.0     wi 23.9.97 bug fixed
%      1.1     wi 26.9.97 latest fit on experimental data was added (option 7)
%      1.2     wi 13.10.97 option 8 added, adapted scattering of a shell/sphere to note 9/ver2 
%      1.3     wi  4.11.97 option 9, 10 and 11 added 
%      1.4     wi 27.05.98 born approximation added (borna.m)
%                 12.03.07 bug in iborn shperes fixed (rtt)
%
%   Uses:
%       borna, ro2epsd, mixmod
%
%
%
%   Copyright (c) 1997 by the Institute of Applied Physics,
%   University of Bern, Switzerland

% constants
c = 2.99; 
roair = 0.001293;
roice = 0.917;
% specular component of scattering coefficient
% usually 0 can be important in new snow!
dgb0h = 0;
dgb0v = 0;
% aus der Theorie scattering coefficient
k = freq*(2*pi/0.299793);
eice = 3.18;
vfi = roi./roice;

% choose the scattering algorithm that should be used
wahl = sccho;


[epsi,epsii] = ro2epsd(roi,Ti,freq);
[epsi,epsii] = mixmod(freq,Ti,Wi,epsi,epsii);


% 6-flux scattering coefficient
if wahl == 1
   gs6 = ((130 * ((freq/50)^2.7)) .* pci.^3) ./ (roi.^1.3 + 0.001);
end


%fit vom 26.8.97 auf alle Daten v-pol, > 11 GHz
if wahl == 2
   gs6 = 0.0704 * (pci.^2.32).*(freq.^1.68).*roi.^(-0.63);
end

% for spheres: Mätzler, J. Appl. Phys. 83(11) 6111-6117 eqs 27+32 iborn
epseff = (2-eice+3.*vfi.*(eice-1)+ sqrt((2-eice+3.*vfi.*(eice-1)).^2+8.*eice))./4;
sphe = (3/32).*(0.001.*pci).^3.*k.^4.*vfi.*(1-vfi).*abs((2.*epseff+1).*(eice-1)./(2.*epseff+eice)).^2;
if wahl == 4
   gs6 = sphe;
end

% for shells(new and recrystalized snow): 
% Mätzler, J. Appl. Phys. 83(11) 6111-6117 eq 39 iborn
epseff = 1+(vfi.*(eice-1).*(2+1/eice))./(3-vfi.*(1-1/eice));
shel = abs(2/3 + 1./(3.*eice.^2)).*(0.001.*pci).*k.^2.*vfi.*(1-vfi).*(eice-1).^2./(16.*epseff);
if wahl == 5
   gs6 = shel;
end

% as linearcombination
if wahl == 6
   a = 0.1664;
   b = 0.2545;
   gs6 = a.*sphe+b.*shel;
end

%fit vom 26.9.97
if wahl == 7
   gs6 = 73.21 * (pci.^3).*((freq./50).^2.68).*roi.^(-1);
end
%fit vom 13.10.97
if wahl == 8
   gs6 = 136 .* (pci.^2.85) .* ((freq./50).^2.5) ./ (roi + 0.001);
end

%fit vom 4.11.97 (without density)
if wahl == 9
   gs6 = 564 .* (pci.^3.0) .* ((freq./50).^2.5);
end

%fit vom 4.11.97 (without density, uses corr. length from exp. fit!)
if wahl == 10
   gs6 = (3.16 .* pci + 295 .* (pci.^2.5)).* ((freq./50).^2.5);
end

%fit vom 4.11.97 (with density, uses corr. length from exp. fit!)
if wahl == 11
   gs6 = (9.20 .* pci - 1.23 .* roi + 0.54).^2.5 .* ((freq./50).^2.5);
end

omega = sqrt((epsi - 1)./epsi);


%Born Approximation
if wahl == 12
  kp = bornsnk(roi,0);
  [gb6,gc6,gf6,gs6] = borna(k,vfi,pci,epsi,eice,epseff,kp);
else
  gb6 = 0.5 .* gs6 .* (1-omega);
  gc6 = 0.25 .* gs6 .* omega;
end

gbiv = zeros(size(gs6));
gbih = zeros(size(gs6));
gtr = zeros(size(gs6));
ga2 = zeros(size(gai));
%dgb2h = zeros(size(gs6));
%dgb2h = dgb2h + 0.25;
%dgb2v = zeros(size(gs6));
%dgb2v = dgb2v + 0.1;

% -> 2 Flux
gtr = (4 .* gc6) ./ (gai + 2 .* gc6);
ga2i = gai .* (1 + gtr);

gbih = (gb6 + dgb0h) + gtr .* gc6;
gbiv = (gb6 + dgb0v) + gtr .* gc6;
