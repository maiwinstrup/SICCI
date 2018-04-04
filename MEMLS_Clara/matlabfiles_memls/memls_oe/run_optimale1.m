function run_optimale1
debug_on_warning(0);
warning('off','all');

data=load('testdata20100419.txt');
%data=load('testdata_marapr.txt');
n=rows(data);
for i=1:39
    tbrow=data(i,:);
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
    [pv,sdia,dtb,pvlc,dtblc]=optimale([Tb]);

    GR=(Tb37v-Tb18v)./(Tb37v+Tb18v);
    sdmc=0.01.*(2.9-782.0.*GR);
    SDlise=1.7701 + 0.017462 .* Tb6v - 0.02801 .* Tb18v + 0.0040926 .*Tb37v;
    %odata=[pv(1) pv(2) pv(3) pv(4) pv(5) dtb(1) dtb(2) dtb(3) dtb(4) dtb(5) dtb(6) dtb(7) dtb(8) dtb(9) dtb(10)]
    %odata=[tbrow(3) tbrow(4) tbrow(5) tbrow(6) pv(1) pv(2) pv(3) pv(4) pv(5) sdmc dtb(1) dtb(2) dtb(4) dtb(5) dtb(9) dtb(10)];
    cost=sqrt(sum(dtb.^2));
    lowcost=sqrt(sum(dtblc.^2));
    %odata=[tbrow(1) tbrow(2) tbrow(3) tbrow(4) tbrow(5) tbrow(6) pv(1) pv(2) pv(3) pv(4) pv(5) sdmc SDlise dtb(1) dtb(2) dtb(3) dtb(4) dtb(5) dtb(6) dtb(7) dtb(8) dtb(9) dtb(10) cost pvlc(1) pvlc(2) pvlc(3) pvlc(4) lowcost];
    odata=[tbrow(3) pv(4) sdmc SDlise lowcost];
    dlmwrite('out_m2.txt',transpose(odata(:)),'-append','precision','%.2f')
    disp(odata(1:5));
end %for

end
