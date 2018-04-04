function [pcc]=pccdistr2(hs,nsnow,nice)
%compute the scatter size
snow_top=0.07; %the correlation length of new snow [mm]
snow_bottom=0.3; %the correlation length of old snow [mm]
icepcc=1.5; %the correlation length of multiyear ice [mm]

%denne beregning aendrer spredningen for tynd sne, den burde foregå et andet sted

ds=(snow_bottom-snow_top)./nsnow;

nelements=nsnow+nice;
snow_layer_thickness=hs./nsnow;
pcc=ones(1,nelements);
pcc(1)=snow_top;
   for i=2:nsnow
      mid_layer_depth=snow_layer_thickness.*i-0.5.*snow_layer_thickness;
      pcc(i)=snow_top+3.0.*((hs./0.5)).*mid_layer_depth.*(sqrt(mid_layer_depth)-mid_layer_depth);
   end %endfor

   for i=1:nice
      pcc(i+nsnow)=icepcc;
   end %endfor
%plot(pcc)
end
