function ro = density(temp)

% unit kg/m^3

temp.cal.C = temp.cal.K -273.15;

if temp.cal.K > 425  % h2o2 boiling point

    temp.liquid.h2o2.K = 425;
else
    temp.liquid.h2o2.K = temp.cal.K;
end

ro.h2o2 = 1597 + 0.0784*temp.liquid.h2o2.K - 0.00197*temp.liquid.h2o2.K^2;

A = 0.9998396;
B = 18.224944*0.001;
C = -7.9221*0.000001;
D = -55.44846*10^(-9);
E = 149.7562*10^(-12);
F = -393.2952*10^(-15);
G = 18.15972510^(-3);

if temp.cal.C >= 100

    temp.liquid.h2o.C = 100;
else
    temp.liquid.h2o.C = temp.cal.C;
end    

ro.h2o = (10^3)*(A + B*temp.liquid.h2o.C +C*temp.liquid.h2o.C^2 + D*temp.liquid.h2o.C^3 + E*temp.liquid.h2o.C^4 + F*temp.liquid.h2o.C^5)/(1+G*temp.liquid.h2o.C);

end