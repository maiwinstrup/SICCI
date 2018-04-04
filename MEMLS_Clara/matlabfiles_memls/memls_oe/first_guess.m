function [temperature, salinity, snow_density, snow_depth, ice_thickness]=first_guess(Tb)

GR=(Tb(7)-Tb(5))./(Tb(7)+Tb(5));

%snow_depth=0.01.*(2.9-782.0.*GR);
snow_depth = 1.7701 + 0.017462 .* Tb(1) - 0.02801 .* Tb(5) + 0.0040926 .* Tb(7);
%if (GR < -0.02)
%   snow_depth=0.1854+0.01*(2.9-260.0*GR);
%end %endif
%snow_depth=0.20;
%temperature=2.93*Tb(1)-1.46*Tb(3)+0.58*Tb(10)-239.6;
%temperature=273.15+1.5*((Tb(1)/0.98)-273.15);
temperature = 1.53*Tb(1) - 136.6;
if (temperature > 272.0)
   temperature = 272.0;
end
if (temperature < 220.0)
   temperature = 220.0;
end
salinity=8.5+100.0*GR;
if (salinity <= -0.04)
   salinity=0.0;
end %endif
if (GR > -0.04)
   salinity=7.0;
end  %endif
snow_density=320.0;
ice_thickness=1.0-25.*GR;
if (ice_thickness > 3.5)
   ice_thickness = 3.5;
end %endif
if (ice_thickness < 1.0)
   ice_thickness= 1.0;
end %endif
%ice_thickness=1.5;
end
