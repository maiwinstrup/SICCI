function read_ice_hdf(file)
debug_on_warning(0);
warning('off','all');


data=load('-hdf5',file);
Tb6v=0.1.*double(data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_06V_DAY);
Tb6h=0.1.*double(data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_06H_DAY);

Tb10v=0.1.*double(data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_10V_DAY);
Tb10h=0.1.*double(data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_10H_DAY);

Tb18v=0.1.*double(data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_18V_DAY);
Tb18h=0.1.*double(data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_18H_DAY);

Tb23v=0.1.*double(data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_23V_DAY);
Tb23h=0.1.*double(data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_23H_DAY);

Tb36v=0.1.*double(data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_36V_DAY);
Tb36h=0.1.*double(data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_36H_DAY);

Tb89v=0.1.*double(data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_89V_DAY);
Tb89h=0.1.*double(data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_89H_DAY);

icecon=0.01.*double(data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_ICECON_DAY);
icediff=data.HDFEOS.GRIDS.NpPolarGrid25km.Data_Fields.SI_25km_NH_ICEDIFF_DAY;

    idx=zeros(304,448);
    land=zeros(304,448);
    spv=zeros(304,448);
    ssd=zeros(304,448);
    spvlc=zeros(304,448);

    %for i=50:250
    for i=1:304
       %for j=50:400
       for j=1:448
          if (Tb6v(i,j) < Tb6h(i,j) || Tb10v(i,j) < Tb10h(i,j) || Tb18v(i,j) < Tb18h(i,j) || Tb36v(i,j) < Tb36h(i,j) || Tb89v(i,j) < Tb89h(i,j))
              continue
          endif
          if (icecon(i,j) < 1.05 && icecon(i,j) > 0.95)
              Tb=[Tb6v(i,j) Tb6h(i,j) Tb10v(i,j) Tb10h(i,j) Tb18v(i,j) Tb18h(i,j) Tb36v(i,j) Tb36h(i,j) Tb89v(i,j) Tb89h(i,j)];
              [pv,sdia,dtb,pvlc,dtblc]=optimaler([Tb]);
              idx(i,j)=1.0;
              spv(i,j)=pv(2);
              ssd(i,j)=sdia(2);
          endif
          if (icecon(i,j) < 1.22 && icecon(i,j) > 1.18)
              land(i,j)=1.0;
          endif
       endfor
    endfor

    GR=(Tb36v-Tb18v)./(Tb36v+Tb18v);
    sdmc=0.01.*(2.9-782.0.*GR);
    GR0610v=(Tb10v-Tb6v)./(Tb10v+Tb6v);
    GR0618v=(Tb18v-Tb6v)./(Tb18v+Tb6v);
    GR0618h=(Tb18h-Tb6h)./(Tb18h+Tb6h);
    GR0636v=(Tb36v-Tb6v)./(Tb36v+Tb6v);
    SDsim= 1.0.*(-0.37735189 + 3.8681547.*GR0610v - 6.7846170.*GR0618v + 20.348159.*GR0618h + 0.27384655.*GR0636v + 0.045854709.*Tb6h - 0.045823858.*Tb18h + 0.0017964592.*Tb36v);

whos
climits=[0.0,0.010];
%imagesc((sdmc.*idx)',climits);
%h=imagesc((spv+land)',climits);
h=imagesc((ssd+land)',climits);
h=colorbar("east");
saveas(h,'reg13_snowd_20170419.png');
%hist(spv,10);
end
