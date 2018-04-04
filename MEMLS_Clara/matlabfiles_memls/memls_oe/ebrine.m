function [eb,ebi]=ebrine(freq,eb0,relax,sig)

% permittivity and loss of brine Ulaby et al. 1986 vol. III E64a
% freq: em frequency
% eb0: static dielectric constant of brine
% relax: relaxation time of brine
% sig: conductivity of brine

f=freq .*1.0e9;
epsiwoo=4.9;
e0=8.854e-12;

eb=epsiwoo+((eb0-epsiwoo) ./(1.+((relax .*f) .^2)));
eb=max(0.0,eb);
ebi=((relax .*f .*(eb0-epsiwoo)) ./(1.+((relax .*f) .^2)))+(sig ./(6.28 .*f .*e0));
ebi=max(0.0,ebi);

