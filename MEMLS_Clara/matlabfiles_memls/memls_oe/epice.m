function epui=epice(T,freq)

% Dielectric constant of pure ice: Mätzler, 1998, Microwave properties of ice and snow, in:
% B. Smitht et al. (eds.) solar system ices, 241-257, Kluwer.

% T: thermometric temperature in K
% freq: Frequency in GHz

epi=3.1884+9.1e-4 .*(T-273);

% The Hufford model for the imagenary part.

theta=300 ./T;
alpha=(0.00504+0.0062 .*theta) .*exp(-22.1 .*theta);
beta=(0.502-0.131 .*theta ./ (1+theta)) .*1e-4 +(0.542e-6 .*((1+theta) ./(theta+0.0073)).^2);

epii=(alpha ./freq) + (beta .*freq);

epui=epi+epii*i;
