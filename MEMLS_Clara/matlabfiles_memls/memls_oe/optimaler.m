function [P_final, S_diag, dTb, Plowcost, dtblowcost] = optimaler(Tb)
    %input parametre is :
    %maalte straalingstemperature Tb [K], 1xn vektor
    %kovariansmatricen Se for Tb nxn matrice
    %foerste gaet paa fysiske vaerdier P0, 1xm vektor
    %kovarians matricen Sp for de fysiske vaerdier, mxm matrice

    %todo: checks on the matrix operations to check for division by zero etc.
    %for every exception raise a specific flag
    %check or print the Tb-Ta
    %check pinv instead of inv

             %read in the Tb's ex. data=load('Tb.txt')
             %Tb=[Tb6v   Tb6h   Tb10v  Tb10h  Tb18v  Tb18h  Tb36v  Tb36h  Tb89v  Tb89h]
             %Tb=[246.44 227.68 242.64 224.19 234.64 220.45 205.37 193.45 169.70 162.65];

    %kovariansmatricen for de estimerede straalingstemperatur vaerdier S
              Se=[0.2 0 0 0 0 0 0 0 0 0;
                  0 0.2 0 0 0 0 0 0 0 0;
                  0 0 1.2 0 0 0 0 0 0 0;
                  0 0 0 1.2 0 0 0 0 0 0;
                  0 0 0 0 1.2 0 0 0 0 0;
                  0 0 0 0 0 1.2 0 0 0 0;
                  0 0 0 0 0 0 1.2 0 0 0;
                  0 0 0 0 0 0 0 1.2 0 0;
                  0 0 0 0 0 0 0 0 20.2 0;
                  0 0 0 0 0 0 0 0 0 20.2];

    %kovariansmatricen for de estimerede fysiske vardier P
               %Sp=[23.1 0.0 0.0 0.0 0.0;
               %    0.0 5.0 0.0 0.0 0.0;
               %    0.0 0.0 15.0 0.0 0.0;
               %    0.0 0.0 0.0 0.01 0.0;
               %    0.0 0.0 0.0 0.0 0.24];
               Sp=[20.0 0.0 0.0;
                   0.0 0.2 0.0;
                   0.0 0.0 0.5];
               %Sp=[24.2 -0.1 -0.1;
               %    -0.1 0.3 0.1;
               %    -0.1 0.1 0.5];
               %Sp=[0.2 0.0 0.0 0.0 0.0;
               %    0.0 0.2 0.0 0.0 0.0;
               %    0.0 0.0 0.2 0.0 0.0;
               %    0.0 0.0 0.0 0.2 0.0;
               %    0.0 0.0 0.0 0.0 0.0];

    %if (det(Sp) == 0 || det(Se) == 0)
    %   continue
    %end

    %vector med 3 fysiske vaerdier
    SDlise=1.7701 + 0.017462.*Tb(1) - 0.02801.*Tb(5) + 0.0040926.*Tb(7);
    snow_depth=0.11;
    snow_depth=SDlise;
    ice_thickness=2.3;
    temperature=230.0;
    P0 = [temperature snow_depth ice_thickness]; %initial vaerdier for isen 

    %oevre og nedre graenser til de fysiske parametre
    L=[190.0 273.15 0.01 1.0 0.5 4.00];
    len_p = length(P0);
    len_tb = length(Tb);
    melements = [len_p len_tb];
    M = ones(melements);
    n=1;
    %indsat med konvergenskriteria
    costmi=100.0;
    cost_threshold=25.0;
    %start while loop
    while n<11 %the simple criteria
    %another criteria
    %while (cost >= cost_threshold)
        n=n+1;
        %kald forward model med seneste gaet. samtidig reference til adjoint model.
        if n==2
            Ta=run_reg(P0(1),P0(2),P0(3));
        else
            %the simulated Tb 
            Ta=run_reg(P(1),P(2),P(3));
        end %if
        dtb=Tb-Ta;
        %start for loop over parametre
        for i=1:len_p
            %beregn partiell afledede ->kald forwardmodel med pertubationer for hver parameter, pertubationer 1%
            if n==2
                %M matricen "the adjoint" foeste iteration bruger initialvaerdierne
                M(i,:) = (run_reg((0.01.*(i==1).*P0(1))+P0(1),(0.01.*(i==2).*P0(2))+P0(2),(0.01.*(i==3).*P0(3))+P0(3)) - ...
                         (run_reg(P0(1),P0(2),P0(3)))) / ...
                         ((0.01*(i==1)*P0(1))+(0.01*(i==2)*P0(2))+(0.01*(i==3)*P0(3)));
                %efterfoelgende iterationer
            else
                M(i,:) = (run_reg((0.01.*(i==1).*P(1))+P(1),(0.01.*(i==2).*P(2))+P(2),(0.01.*(i==3).*P(3))+P(3)) - ...
                         (run_reg(P(1),P(2),P(3)))) / ...
                         ((0.01.*(i==1).*P(1))+(0.01.*(i==2).*P(2))+(0.01.*(i==3).*P(3)));
            end %if
        end %for
        M = transpose(M);
        %beregn S
        %if (det(inv(Sp) + transpose(M)*inv(Se)*M) == 0.0)
        %   continue
        %end
        S = inv(inv(Sp) + transpose(M)*inv(Se)*M);
        %check konvergens kriteria f.eks. n>5
        %hvis n>5 (antallet af iterationer) print S diagonal + alt andet
        %if (n>20 | cost < cost_threshold)
        if (n>10)
            P_final = P;
            S_diag = diag(S);
            dTb=Tb-Ta;
            %indsaet sammen med konvergenskriteria
            break
        elseif n==2
            %transpose frem og tilbage, er det nodvendigt?
            P0 = transpose(P0);
            Tb = transpose(Tb);
            Ta = transpose(Ta);
            %estimer fysiske vaerdier ud fra foerste gaet
            %print statement
            P = P0 + S*(transpose(M)*inv(Se)*(Tb - Ta) + inv(Sp)*(P0 - P0));
            %Det virker strengt taget ikke til at vaere noedvendigt med fysiske begraensninger, for aabent vand
            if (P(1) < L(1)) 
                P(1) = L(1)+0.1;
            end %if
            if (P(1) > L(2)) 
                P(1) = L(2);
            end %if
            if (P(2) < L(3)) 
                P(2) = L(3);
            end %if
            if (P(2) > L(4)) 
                P(2) = L(4);
            end %if
            if (P(3) < L(5)) 
                P(3) = L(5)+0.1;
            end %if
            if (P(3) > L(6)) 
                P(3) = L(6);
            end %if

            P_1 = P;
            %transpose frem og tilbage, er det nodvendigt?
            M = transpose(M);
            P0 = transpose(P0);
            Tb = transpose(Tb);
            Ta = transpose(Ta);
            P = transpose(P);
            dtb=Tb-Ta;
            cost=sqrt(sum(dtb.^2));
        else
            %transpose frem og tilbage, er det nodvendigt?
            P0 = transpose(P0);
            Tb = transpose(Tb);
            Ta = transpose(Ta);

            %2. 3. 4. iteration frem mod fysiske vaerdier der faar Tb til at passe
            %aendre til p0
            %P = P0 + S*(transpose(M)*inv(Se)*(Tb - Ta) + inv(Sp)*(P0 - P_1));
            P = P_1 + S*(transpose(M)*inv(Se)*(Tb - Ta) + inv(Sp)*(P0 - P_1));
            if (P(1) < L(1)) 
                P(1) = L(1)+0.1;
            end
            if (P(1) > L(2)) 
                P(1) = L(2);
            end
            if (P(2) < L(3)) 
                P(2) = L(3);
            end %if
            if (P(2) > L(4)) 
                P(2) = L(4);
            end %if
            if (P(3) < L(5)) 
                P(3) = L(5)+0.1;
            end %if
            if (P(3) > L(6)) 
                P(3) = L(6);
            end %if

            P_1 = P;
            %transpose frem og tilbage, er det nodvendigt?
            M = transpose(M);
            P0 = transpose(P0);
            Tb = transpose(Tb);
            Ta = transpose(Ta);
            P = transpose(P);
            dtb=Tb-Ta;
            cost=sqrt(sum((Tb-Ta).^2));
            if cost < costmi
               costmi = cost;
               Plowcost = P;
               dtblowcost = dtb; 
            %tilbage til while loop start
        end %if
    end %while
end
