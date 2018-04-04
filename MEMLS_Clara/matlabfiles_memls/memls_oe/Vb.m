function volbrine=Vb(T,sal)

% the volume of brine in sea ice Ulaby et al. 1986, E71
% T: thermometric temperature [C]
% sal: sea ice salinity [ppt]
if (abs(T) < 0.1) 
    volbrine = 0.001 .*sal .*49.717;
    return
end %endif
volbrine = 0.001 .*sal .*((-49.185 ./T)+0.532);
end