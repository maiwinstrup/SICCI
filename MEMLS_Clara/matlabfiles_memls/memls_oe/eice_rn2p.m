function emi=eice_rn2p(eh,ei,vi)

% Polder / Van Santen dielectric constant of sea ice with random brine needles
% e.g. Shokr (1998) Arctic sea ice in the microwave c-band, IEEE TGRS 36(2), 463-78.

measureA=1.2;
measureB=1.2;

% initial estimate: Vant et al. (1978), J. Appl. Phys. 49(3), 1264-80.
emi=(3.05+0.72 .*vi)+(0.024+0.33 .*vi) .*i;

m=0;

% iterate to solve for emi: the dielectric const. of the mixture
while measureA > 0.001 || measureB > 0.001

 f1=(ei-eh) ./ (ei+emi);
 f2=5. .*emi+ei;
 est=eh+((vi ./ 3.) .* f1 .* f2);

 measureA=abs(real(emi)-real(est));
 measureB=abs(imag(emi)-imag(est));

 emi=est;
 m=m+1;

 if (m>30) break; end
end

emi;
