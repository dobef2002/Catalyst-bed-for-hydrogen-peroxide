function mu = viscosity(temp,R,X,ndot,mw)

temp.cal.C =temp.cal.K - 273.15;

massfra.liquid.h2o2 = (ndot.h2o2*mw.h2o2)/(ndot.h2o2*mw.h2o2+ndot.h2o*mw.h2o);
massfra.liquid.h2o = (ndot.h2o*mw.h2o)/(ndot.h2o2*mw.h2o2+ndot.h2o*mw.h2o);

if temp.cal.C < 150
   mu.liquid.h2o2 = 0.001 * 0.2544152386* exp(332762538.8/(R*temp.cal.K^3)); % cpoise

else
   mu.liquid.h2o2 = 0.001 * 0.2544152386* exp(332762538.8/(R*(150+273.15)^3)); % cpoise
end

if temp.cal.C < 100
   B = (1.3272*(20-temp.cal.C) - 0.001053*(temp.cal.C -20)^2)/(temp.cal.C + 105); 
   mu.liquid.h2o = 1.002*0.001*10^B;

else
   B = (1.3272*(20-100) - 0.001053*(100 -20)^2)/(100 + 105); 
   mu.liquid.h2o = 1.002*0.001*10^B;
end

mu.liquid.mixture = mu.liquid.h2o2*(massfra.liquid.h2o2)^2 + 2*massfra.liquid.h2o2*massfra.liquid.h2o*(mu.liquid.h2o2*mu.liquid.h2o)^0.5 + mu.liquid.h2o*(massfra.liquid.h2o)^2;

mu.gas.h2o2 = 10^(-7)*(134 + 0.35*(temp.cal.C - 100)-14);
mu.gas.h2o = 10^(-7)*((39.37*temp.cal.K^1.5)/(3315 - temp.cal.K + 0.001158*tmp.cal.K^2));

mu.gas.o2 = 10^(-3)*0.02018*((0.555*526.05+127)/(0.555*(temp.cal.K*1.8)+127))*((temp.cal.K*1.8/526.05)^1.5);
mu.gas.mixture = (mu.gas.h2o2*X.gas.h2o2*(mw.h2o2)^0.5 + mu.gas.h2o*X.gas.h2o*(mw.h2o)^0.5 + mu.gas.o2*X.gas.o2*(mw.o2)^0.5) / (X.gas.h2o2*(mw.h2o2)^0.5 + X.gas.h2o*(mw.h2o)^0.5 + X.gas.o2*(mw.o2)^0.5);



end


