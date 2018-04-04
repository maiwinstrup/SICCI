function run_optimaler1
debug_on_warning(0);
warning('off','all');

data=load('testdata20100419.txt');
%data=load('testdata_marapr.txt');
n=rows(data);
%39?
for i=1:39
    tbrow=data(i,:);
    snow=tbrow([3]);
    ice=tbrow([5]);
    Tb=tbrow([7:16]);
    Tb6v=Tb(1);
    Tb6h=Tb(2);
    Tb10v=Tb(3);
    Tb10h=Tb(4);
    Tb18v=Tb(5);
    Tb18h=Tb(6);
    Tb37v=Tb(7);
    Tb37h=Tb(8);
    Tb89v=Tb(9);
    Tb89h=Tb(10);
    if (Tb6v < Tb6h || Tb10v < Tb10h || Tb18v < Tb18h || Tb37v < Tb37h || Tb89v < Tb89h)
       continue
    end
    for j=1:10
       if (Tb(j) < 100.0 || Tb(j) > 273.15)
          continue
       end %if
    end %for
    [pv,sdia,dtb,pvlc,dtblc]=optimaler([Tb]);

    GR=(Tb37v-Tb18v)./(Tb37v+Tb18v);
    sdmc=0.01.*(2.9-782.0.*GR);
    GR0610v=(Tb10v-Tb6v)./(Tb10v+Tb6v);
    GR0618v=(Tb18v-Tb6v)./(Tb18v+Tb6v);
    GR0618h=(Tb18h-Tb6h)./(Tb18h+Tb6h);
    GR0636v=(Tb37v-Tb6v)./(Tb37v+Tb6v);
    %Lise's equations
    SDsim_lise = 1.7701 + 0.017462.*Tb6v - 0.02801.*Tb18v + 0.0040926*Tb37v;
    Tsi6_lise = 1.144.*Tb6v - 0.815/SDsim_lise - 27.08;
    Tsi10_lise = 1.101.*Tb10v - 0.999/SDsim_lise - 14.28;
    Tsi_lise = 0.5.*(Tsi6_lise + Tsi10_lise);
    ks=0.3;
    ki=2.1;
    hi=ice;
    Ta=Tb89v+15;
    Tw=271.35;
    hs=(hi.*(ks.*Tsi_lise-ks.*Ta))./(ki.*Tw-ki.*Tsi_lise);

    SDsim= -0.377351898662 + 3.86815473658.*GR0610v -6.78461704699.*GR0618v+ 20.3481598581.*GR0618h+ 0.273846558015.*GR0636v+ 0.0458547097079.*Tb6h -0.0458238581105.*Tb18h+ 0.00179645929531.*Tb37v;
    SDlise=1.7701 + 0.017462.*Tb6v - 0.02801.*Tb18v + 0.0040926.*Tb37v;
    %odata=[pv(1) pv(2) pv(3) dtb(1) dtb(2) dtb(3)]
    %odata=[tbrow(3) tbrow(4) tbrow(5) tbrow(6) pv(1) pv(2) pv(3) pv(4) pv(5) sdmc dtb(1) dtb(2) dtb(4) dtb(5) dtb(9) dtb(10)];
    cost=sqrt(sum(dtb.^2));
    lowcost=sqrt(sum(dtblc.^2));
%    odata=[tbrow(1) tbrow(2) tbrow(3) tbrow(4) tbrow(5) tbrow(6) pv(1) pv(2) pv(3) sdmc dtb(1) dtb(2) dtb(3) dtb(4) dtb(5) dtb(6) dtb(7) dtb(8) dtb(9) dtb(10) cost pvlc(1) pvlc(2) pvlc(3) lowcost SDsim];
    %odata=[tbrow(1) tbrow(2) tbrow(3) tbrow(4) tbrow(5) tbrow(6) pv(1) pv(2) pv(3) sdmc cost pvlc(1) pvlc(2) pvlc(3) lowcost SDsim];
    odata=[snow ice hs sdmc cost lowcost SDsim pv(2) SDlise];
    dlmwrite('out_testdata.txt',transpose(odata(:)),'-append','precision','%.2f')
    disp(odata(1:9));
end %for

end
