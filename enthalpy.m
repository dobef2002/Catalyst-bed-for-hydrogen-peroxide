function h = enthalpyT(temp,h,T)  % calculate enthalpy (kJ/kmol)

h.hf.liquid.h2o2 = -187860;
h.hf.gas.h2o2 = -136310;
h.hf.liquid.h2o = -285830;
h.hf.gas.h2o = -241820;
h.hf.gas.o2 = 0;

temp.liquid.h2o2.K = T;
h.dhs.liquid.h2o2 = 89.377*(temp.liquid.h2o2.K - temp.ref);

h.enthalpy.liquid.h2o2 = h.hf.liquid.h2o2 + h.dhs.liquid.h2o2;

const.A.gas.h2o2 = 34.25667;
const.B.gas.h2o2 = 55.18445;
const.C.gas.h2o2 = -35.15443;
const.D.gas.h2o2 = 9.087440;
const.E.gas.h2o2 = -0.422157;
const.F.gas.h2o2 = -149.9098;
const.G.gas.h2o2 = -136.1064;

const.A.liquid.h2o = -203.606;
const.B.liquid.h2o = 1523.26;
const.C.liquid.h2o = -3196.413;
const.D.liquid.h2o = 2474.455;
const.E.liquid.h2o = 3.855326;
const.F.liquid.h2o = -256.54786;
const.G.liquid.h2o = -285.8304;

const.A.gas.h2o = 30.09200;
const.B.gas.h2o = 6.832515;
const.C.gas.h2o = 6.7934535;
const.D.gas.h2o = -2.534480;
const.E.gas.h2o = 0.082139;
const.F.gas.h2o = -250.8810;
const.G.gas.h2o = -241.8264;

const.A.gas.o2.Tle700 = 31.32234;
const.B.gas.o2.Tle700 = -20.23531;
const.C.gas.o2.Tle700 = 57.86644;
const.D.gas.o2.Tle700 = -36.50624;
const.E.gas.o2.Tle700 = -0.007374;
const.F.gas.o2.Tle700 = -8.903471;
const.G.gas.o2.Tle700 = 0;

const.A.gas.o2.Tge700 = 30.03235;
const.B.gas.o2.Tge700 = 8.772972;
const.C.gas.o2.Tge700 = -3.988133;
const.D.gas.o2.Tge700 = 0.788313;
const.E.gas.o2.Tge700 = -0.741599;
const.F.gas.o2.Tge700 = -11.32468;
const.G.gas.o2.Tge700 = 0;

temp.liquid.h2o.K = T;
t = T/1000;
t_L_h2o = temp.liquid.h2o.K/1000;

h.dhs.gas.h2o2 = 1000*(const.A.gas.h2o2*t + 0.5*const.B.gas.h2o2*t^2 + (1/3)*const.C.gas.h2o2*t^3 + 0.25*const.D.gas.h2o2*t^4 - (const.E.gas.h2o2/t) + const.F.gas.h2o2 - const.G.gas.h2o2);
h.dhs.liquid.h2o = 1000*(const.A.liquid.h2o*t_L_h2o + 0.5*const.B.liquid.h2o*t_L_h2o^2 + (1/3)*const.C.liquid.h2o*t_L_h2o^3 + 0.25*const.D.liquid.h2o*t_L_h2o^4 - (const.E.liquid.h2o/t_L_h2o) + const.F.liquid.h2o - const.G.liquid.h2o);
h.dhs.gas.h2o = 1000*(const.A.gas.h2o*t + 0.5*const.B.gas.h2o*t^2 + (1/3)*const.C.gas.h2o*t^3 + 0.25*const.D.gas.h2o*t^4 - (const.E.gas.h2o/t) + const.F.gas.h2o - const.G.gas.h2o);

if temp.cal.K <= 700

h.dhs.gas.o2 = 1000*(const.A.gas.o2.Tle700*t + 0.5*const.B.gas.o2.Tle700*t^2 + (1/3)*const.C.gas.o2.Tle700*t^3 + 0.25*const.D.gas.o2.Tle700*t^4 - (const.E.gas.o2.Tle700/t) + const.F.gas.o2.Tle700 - const.G.gas.o2.Tle700);

else

h.dhs.gas.o2 = 1000*(const.A.gas.o2.Tge700*t + 0.5*const.B.gas.o2.Tge700*t^2 + (1/3)*const.C.gas.o2.Tge700*t^3 + 0.25*const.D.gas.o2.Tge700*t^4 - (const.E.gas.o2.Tge700/t) + const.F.gas.o2.Tge700 - const.G.gas.o2.Tge700);

end



h.enthalpy.gas.h2o2 = h.hf.gas.h2o2 + h.dhs.gas.h2o2;
h.enthalpy.liquid.h2o = h.hf.liquid.h2o + h.dhs.liquid.h2o;
h.enthalpy.gas.h2o = h.hf.gas.h2o + h.dhs.gas.h2o;
h.enthalpy.gas.o2 = h.hf.gas.o2 + h.dhs.gas.o2;

end

