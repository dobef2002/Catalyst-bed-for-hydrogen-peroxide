function pressure = p_saturationT(temp,pressure,T)

temp.liquid.h2o2.K = T;

A = -4025.3;
B = 12.966;
C = 4.6055e-3;
D = 44.576;

order = A/temp.liquid.h2o2.K + B*log10(temp.liquid.h2o2.K) + C*temp.liquid.h2o2.K + D;

pressure.sat.h2o2 = 133.322365*(10^order);

temp.liquid.h2o.K = T;

A = -2892.3698;
B = -2.892736;
C = -4.9369728e-3;
D = -5.606905e-6;
E = -4.645869e-9;
F = -3.7874e-12;
G = 19.301141;

order = A*temp.liquid.h2o.K + B*log10(temp.liquid.h2o.K) + C*temp.liquid.h2o.K + D*temp.liquid.h2o.K^2 + E*temp.liquid.h2o.K^3 + F*temp.liquid.h2o.K^4 +G;

pressure.sat.h2o = 133.322365*(10^order);

end