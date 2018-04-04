function [data]=memls_mod(num,di,Ti,Wi,roi,pci,sal,type,si)

%   MEMLS
%   main program
%   of a snowpack/and packice by frequency.
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland
%   Extension to sea ice by rtt dmi.dk

% her kunne man godt beregne reflektionskoefficienterne som funktion af temperaturen
teta=53;
s0h=0.75;
s0v=0.25;
Tsky=0;
Tgnd=271.35;
sccho=4; %4: iborn for spheres, 5: iborn for shells and spheres, now a function of type.
teta = (teta * pi) / 180;

N = max(size(num));

if N == 0 
   return
end

roi = roi./1000;

%  loop over frequency
freqs= [1.4, 6.9, 10.7, 18.7, 23.8, 36.5, 50.0, 89.0, 157.0, 183.0, 243.0, 325.0, 448.0, 664.0];
s0hs= [0.40, 0.39, 0.39, 0.38, 0.37, 0.35, 0.32, 0.26, 0.19, 0.18, 0.18, 0.18, 0.18, 0.18];
s0vs= [0.39, 0.37, 0.37, 0.35, 0.33, 0.30, 0.28, 0.21, 0.15, 0.13, 0.13, 0.13, 0.13, 0.13];

x  = freqs;
yh = freqs;
yv = freqs;
yh0 = freqs;
yv0 = freqs;
yh100 = freqs;
yv100 = freqs;
tb = [freqs,freqs]; % To hold the tbs in runalg sequence
emis = [freqs,freqs]; % To hold the tbs in runalg sequence

      for ifreq = 1:14
          freq=freqs(ifreq);
          s0h=s0hs(ifreq);
          s0v=s0vs(ifreq);
          [epsi,epsii] = ro2epsd(roi,Ti,freq);
          [epsi,epsii] = mixmod(freq,Ti,Wi,epsi,epsii);
%   select dielectric scheme for FY or MY ice
          fy=(type==3);
          my=(type==4);
          [epsi,epsii] = sie(fy,sal,Ti,freq,epsi,epsii);
          [epsi,epsii] = mysie(my,roi,Ti,sal,freq,epsi,epsii);
          [gai] = abscoeff(epsi,epsii,Ti,freq,Wi);
          ns = sqrt(epsi);
          tei = [asin(sin(teta)./ns);teta];
          dei = pfadi(tei,di);
          [sih,siv] = fresnelc(tei,[epsi;1]);
          [rnum,rroi,repsi,repsii,rtei,rsih,rsiv,rdi,rdei,rTi,rpci,rWi,rgai,rtype,rsal] = slred(num,roi,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,freq,Wi,gai,type,sal);
          
          [gbih,gbiv,gs6,ga2i] = sccoeff(rroi,rTi,rpci,freq,rWi,rgai,sccho);
	  rmeteo=(rtype==1);
	  rrecry=(rtype==2);
	  [gbih,gbiv,gs6,ga2i] = meteo_sc(rroi,rTi,rpci,freq,rWi,rgai,rmeteo,gbih,gbiv,gs6,ga2i);
	  [gbih,gbiv,gs6,ga2i] = recry_sc(rroi,rTi,rpci,freq,rWi,rgai,rrecry,gbih,gbiv,gs6,ga2i);
%  select MY or FY ice scattering modules
          rfy=(rtype==3);
	  rmy=(rtype==4);
          [gbih,gbiv,gs6,ga2i] =    scice(rfy,gbih,gbiv,gs6,ga2i,rTi,rsal,freq,rpci);
          [gbih,gbiv,gs6,ga2i] = scice_my(rmy,gbih,gbiv,gs6,ga2i,rTi,rroi,freq,rpci,rsal);
%  recompute the 2flux absorption coefficient 
          [gbih,gbiv,gs6,ga2i] = absorp2f(gbih,gbiv,gs6,ga2i,repsi,repsii,rroi,rTi,rpci,freq,rWi,rgai);
          [rdei,rtei,tscat] = pfadc(teta,rdi,repsi,gs6);
          rsih = [s0h;rsih];
          rsiv = [s0v;rsiv];
          [rsih,rsiv] = polmix(tscat,rsih,rsiv);
          rtei.*180./pi;
          [ri,ti]  = rt(ga2i,gbih,rdei); 
          
          Dh   = layer(ri,rsih,ti,rTi,Tgnd,Tsky);
          N = max(size(rroi));
          Tbh = (1-rsih(N+1))*Dh(N) + rsih(N+1)*Tsky;
          yh(ifreq) = Tbh; 
          tb(2*ifreq) = Tbh; 
          [freq,N,Tbh];
          [Dh;Tbh];
          [ri,ti]  = rt(ga2i,gbiv,rdei);
          tscat;
          rsih;
          rsiv;
          
          Dv   = layer(ri,rsiv,ti,rTi,Tgnd,Tsky);
          Tbv = (1-rsiv(N+1))*Dv(N) + rsiv(N+1)*Tsky;
          tb(2*ifreq-1) = Tbv;
          [freq,N,Tbv];
          [Dv;Tbv];

%emissivity
          eDh   = layer(ri,rsih,ti,rTi,Tgnd,0);
          eTbh = (1-rsih(N+1))*eDh(N);
	  yh0(ifreq) = eTbh;

          eDh = layer(ri,rsih,ti,rTi,Tgnd,100);
          eTbh = (1-rsih(N+1))*eDh(N) + rsih(N+1)*100;
          yh100(ifreq) = eTbh;

          eDv   = layer(ri,rsiv,ti,rTi,Tgnd,0);
          eTbv = (1-rsiv(N+1))*eDv(N);
          yv0(ifreq) = eTbv;

          eDv   = layer(ri,rsiv,ti,rTi,Tgnd,100);
          eTbv = (1-rsiv(N+1))*eDv(N) + rsiv(N+1)*100;
          yv100(ifreq) = eTbv;

      end % freq

% calculate emissivities
          rv = (yv100 - yv0)./100;
          rh = (yh100 - yh0)./100;
          ev = 1 - rv;
          eh = 1 - rh;
          Teffv = yv0./ev;
          Teffh = yh0./eh;

          x;eh;ev;yh;yv;yh./eh;yv./ev;

          data=real([tb]);

 end %endfunction
