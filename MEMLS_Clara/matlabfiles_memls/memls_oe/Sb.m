function s = Sb(T)

% Salinity of brine
% Vant et al. 1978 J. Appl. Phys. 49(3),1264-1280 & Ulaby et al. 1986 E63
% T: thermometric temperature of brine [C]

if T<-1.8 && T>=-8.2
  s=1.725-18.756 .*T-0.3964 .*(T.^2);
elseif T<-8.2 && T>=-22.9
  s=57.041-9.929 .*T-0.16204 .*(T.^2)-0.002396 .*(T.^3);
elseif T<-22.9 && T>=-36.8
  s=242.94+1.5299 .*T+0.0429 .*(T.^2);
elseif T<-36.8 && T>=-43.2
  s=508.18+14.535 .*T+0.2018 .*(T.^2);
else
%  s=1.725-18.756 .*T-0.3964 .*(T.^2);
%  disp('error in Sb or T outside valid range [-43.2<T[C]<-2]')
  s=34.2;
end

