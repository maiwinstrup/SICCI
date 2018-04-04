function run_optimale
debug_on_warning(0);
warning('off','all');

data=load('AMSR_SIC1_2008.txt');
n=rows(data);
for i=8:16
    tbrow=data(i,:);
    Tb=tbrow([1:6 9:12]);
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
    [pv,sdia,dtb]=optimale([Tb]);

    %odata=[pv(1) pv(2) pv(3) pv(4) pv(5) dtb(1) dtb(2) dtb(3) dtb(4) dtb(5) dtb(6) dtb(7) dtb(8) dtb(9) dtb(10)]
    odata=[pv(1) pv(2) pv(3) pv(4) pv(5) mean(abs(dtb))]
end %for

end
