function rho=dens(TK,S,porosity)

%Eiken, chap 2 + Cox and Weeks

T=TK-273.15;
%porosity in fraction of one
Vv=porosity;

 a = [-4.732,-22.45,-6.397E-1,-1.074E-2];

if (T>-22.9 && T<-2)
   a = [-4.732,-22.45,-6.397E-1,-1.074E-2];
end %endif
if (T<=-22.9)
   a = [9.899E3,1.309E3,5.527E1,7.160E-1];
end %endif

F1 = a(1) + a(2).*T + a(3).*T.^2 + a(4).*T.^3;

Sb=1.725-18.756.*T-0.3964.*T.^2;

rhoi=0.917-1.403E-4.*T;
rhob=1.0+0.0008.*Sb;
Vb=rhoi.*S./(F1+S.*rhoi-S.*rhob); 

if (T<0.0 && T>-2.0)
   ba=[-4.1221E-2,-18.407,5.8402E-1,2.1454E-1];
   bb=[9.0312E-2,-1.6111E-2,1.2291E-4,1.3603E-4];
   F1 = ba(1) + ba(2).*T + ba(3).*T.^2 + ba(4).*T.^3;
   F2 = bb(1) + bb(2).*T + bb(3).*T.^2 + bb(4).*T.^3;
   Vb=rhoi.*S./(F1-rhoi.*S.*F2);
end %endif

rhop=(rhoi.*(1.0-Vb)+rhob.*Vb);
Va=(1.0-(rhop./rhoi))+(rhop.*S./Sb).*(1.0./rhoi - 1.0./rhob);

rho=1000.0.*(rhoi.*(1.0-Vb-Va-Vv)+rhob.*Vb);

end %endfunction
