function [gbih,gbiv,gs6,ga2i] = absorp2f(gbih,gbiv,gs6,ga2i,epsi,epsii,roi,Ti,pci,freq,Wi,gai)

%   calculates the scattering coefficient from structural parameters
%
%   [gbih,gbiv,gs6,ga2i] = absorp2f(gbih,gbiv,gs6,ga2i,roi,Ti,pci,freq,Wi,gai)
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

roi=roi./1000.0;
% constants
c = 2.99; 
roair = 0.001293;
roice = 0.926;
% specular component of scattering coefficient
% usually 0 can be important in new snow!
dgb0h = 0;
dgb0v = 0;
% aus der Theorie scattering coefficient
k = freq*(2*pi/0.299793);
eice = 3.18;
vfi = roi./roice;

omega = sqrt((epsi - 1)./epsi);

  gb6 = 0.5 .* gs6 .* (1-omega);
  gc6 = 0.25 .* gs6 .* omega;

gbiv = zeros(size(gs6));
gbih = zeros(size(gs6));
gtr = zeros(size(gs6));
ga2 = zeros(size(gai));

% -> 2 Flux
gtr = (4 .* gc6) ./ (gai + 2 .* gc6);
ga2i = gai .* (1 + gtr);

gbih = (gb6 + dgb0h) + gtr .* gc6;
gbiv = (gb6 + dgb0v) + gtr .* gc6;
