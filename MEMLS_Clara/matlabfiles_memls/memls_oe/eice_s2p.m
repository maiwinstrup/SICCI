function eeff=eice_s2p(e1,e2,v)

% improved born approximation by C. Mätzler (1998). J. Appl. Phys. 83(11),6111-7
% Polder/VanSanten mixing formulae for spheical inclusions
%effective dielectric constant of medium consisting of e1 and e2
%e1: dielectric constant of background
%e2: dielectric constant of sherical inclusion
%v: fraction of inclusions

eeff=0.25 .*(2 .*e1-e2+3 .*v .*(e2-e1)+sqrt((2 .*e1-e2+3 ...
.*v .*(e2-e1)).^2 +8 .*e1 .*e2));

