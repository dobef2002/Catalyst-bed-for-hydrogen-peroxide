%%% Catalyst bed calculation
%%% solve decompose rate and pressure drop ODE with ODE45

clc;
clear;

%%% USER INPUT

% % propellant properties
R = 8.314; % kJ/(kmol*k) ideal gas constant
mw.h2o2 = 34;
mw.h2o = 18;
mw.o2 = 32;
WT = 0.95; % weight percentage

% % pressure
pessure.inlet = 2500000; % pa
pressure.drop = 0;
pressure.cal = pessure.inlet - pressure.drop;

% % bed dimension
bed.dia = 0.29;             % catalyst bed diameter (m)
bed.length = 0.0294;        % catalyst bed length   (m)
bed.pellet_dia = 0.002;     % catalyst pellet or mesh diameter/wire diameter (m)
bed.porosity = 0.375+0.34*bed.pellet_dia/bed.dia;    % bey equation
bed.segment = 6000;
bed.dx = bed.length/bed.segment;
bed.area = (pi*bed.dia^2)/;
bed.ini = 0;

% % reaction rate properties
reac.Ns= 1; reac.K1 = 1;
reac.Aspecific = =19.7941;                   % Arrhenius pre-exponential factor per unit specific area,[kmol/m2⋅s]
reac.Ef2 =  15000;                           % activation energy, kJ/kmol
reac.as = (1-bed.porosity)6/bed.pellet_dia;  % specific area of catalyst, m^(-1)
reac.Af2 = as*As[ecific;                     % Arrhenius pre-exponential factor, [kmol/m3⋅s]
bed.area_eff = bed.area*(1-bed.porosity)     % effective area

% % mass flow rate kg/s
mdot.total = 0.067;                         %
mdot.h2o2 = mdot.total*WT;
mdor.h2o = mdot.total(1-WT);

% % initial mole flow rate kmol/s

ndot.h2o2 = mdot.h2o2/mw.h2o2;
ndot.h2o = mdot.h2o /mw.h2o;
ndot.o2 = 0;

ndot.initial.h2o2 = ndot.h2o2
ndot.initial.h2o = ndot.h2o;
ndot.initial.o2 = ndot.o2;

% % initial value

temp.enviroment.K = 300; %(k)
temp.enviroment.C = temp.enviroment.K - 273.15; % (C)
temp.initial = temp.enviroment.K;
temp.cal.K = temp.initial;
temp.cal.C = temp.cal.K - 273.15
temp.ref = 298.15;
h = enthalpy(temp);
h.total.initial = ndot.initial.h2o2*h.enthalpy.liquid.h2o2 + ndot.initial.h2o*h.enthalpy.liquid.h2o;


decompose.ndot_initial = 0.0000001*ndot.h2o2;
decompose.ndot = decompose.ndot_initial;
decompose.lamda = 0;

evap_cal = 0;

% % condition check 
checkflag.evap = 0;
checkflag.decomp = 0;
checkflag.div = 0;

% % matrx dimension
pressure_plot = zeros(bed.segment,1);
decompose_ndot_plot = zeros((bed.segment,1);
decompose_lamda_plot = zero(bed.segment,1);
temp_plot = zero(bed.segment,1);
evap_plot = zero(bed.segment,1);

for i = 1:bed.segment

% % mole flow rate kmol/s

ndot.h2o2 = ndot.initial.h2o2 - decompose.ndot;
ndot.h2o = ndot.initial.h2o + decompose.ndot;
ndot.o2 = ndot.initial.o2 + 0.5*decompose.ndot;

if ndot.h2o2 < 0
    disp('H2O2 fully decompose')
    checkflag.decomp = 1;
    ndot.h2o2 = ndot.initial.h2o2 + decompose.ndot;
    ndot.h2o = ndot.initial.h2o - decompose.ndot;
    ndot.o2 = ndot.initial.o2 - 0.5*decompose.ndot;
    decompose.lamda = 100;
elseif (ndot.initial.h2o2 - decompose.ndot) < 0.0001
    disp('H2O2 fully decompose')
    checkflag.decomp = 1;
    decompose.lamda = 100;
end

switch  checkflag.decomp
     case 0

       % % enthalp conservation

          inputs.evap_cal = evap_cal;
          inputs.temp = temp;
          inputs.ndot = ndot;
          inputs.pressure = pressure;
          inputs.h = h;
          inputs.checkflag = checkflag;

                            %  [evap_cal,T_cal,temp] = vpasolveT(inputs);
          outputs = vpasolveT(inputs);
        
          evap_cal = outputs.evap_cal;
          temp = outputs.temp;
          T_cal = outputs.T_cal;

          if evap_cal >= 1
              evap_cal;
              checkflag.evap = 1;
          end

      % % mole fraction
       X.gas.o2 = ndot.o2 / (ndot.o2 +evap_cal * (ndot.h2o2 + ndot.h2o));
       X.gas.h2o2 = ndot.h2o2*evap_cal / (ndot.o2 + evap_cal*(ndot.h2o2 + ndot.h2o));
       X.gas.h2o = ndot.h2o*evap_cal / (ndot.o2 + evap_cal* (ndot.h2o2 + ndot.h2o));
       X.liquid.h2o2 = ndot.h2o2 / (ndot.h2o2 + ndot.h2o);
       X.liquid.h2o = ndot.h2o / (ndot.h2o2 + ndot.h2o);

       % % flow properties calculation

       ro = density(temp);
       volflux.liquid = (1-evap_cal) * ((ndot.h2o2 * mw.h2o2 / ro.h2o2) + (ndot.h2o * mw.h2o / ro.h2o));
       volflux.gas = 0.001 * (ndot.o2 + evap_cal * (ndot.h2o2 + ndot.h2o)) * (R * temp.cal.K / pressure.cal);
       ro.total = (ndot.h2o2 * mw.h2o2 + ndot.h2o * mw.h2o + ndot.o2 * mw.o2) / (volflux.liquid + volflux.gas);
       void = volflux.gas / (volflux.liquid + volflux.gas);

            % % Molar concentration of all species

            concen.h2o2 = ndot.h2o2 / (volflux.liquid + volflux.gas);
            concen.h2o = ndot.h2o / (volflux.liquid + volflux.gas);
            concen.o2 = ndot.o2 / (volflux.liquid + volflux.gas);

            % % Dynamic and kinematic viscosities of the homogeneous flow

            vis_inputs.R = R;
            vis_inputs.temp = temp;
            vis_inputs.ndot = ndot;
            vis_inputs.X = X;
            vis_inputs.mw = mw;
            
            mu = viscosity(vis_inputs);
            mu.mix = mu.liquid.mixture * (1-void) + (mu.gas.mixture * void);  % Dynamic viscosities
            dyv =  mu.mix / ro.total;                                         % kinematic viscosities
            Rep = supV * bed.pellet_dia / dyv;                                % Reynolds number

            % Assume that thee concentration of H2O2 on the catalyst surface is identical to that of the H2O2 in the flow %

            concen.surface.h2o2 = concen.h2o2;

            supV = (volflux.liquid + volflux.gas) / bed.area_eff;
            rdot = Af2 * (exp(-Ef2 / R * temp.cal.K)) * Ns *(K1 * concen.surface.h2o2 / (1 + K1 * concen.surface.h2o2));

            ode_inputs.volflux = volflux;
            ode_inputs.Rep = Rep;
            ode_inputs.supV = supV;
            ode_inputs.rdot = rdot;
            ode_inputs.bed = bed;
            ode_inputs.ro = ro;
                        
            xspan = [bed.ini bed.ini + bed.dx/2 bed.ini + bed.dx];
            A = ode45(@(x,Y) fun_dlamda(x,Y,ode_inputs),xspan,[decompose.ndot,pressure.cal]);
            spot = length(A.x);
            decompose.ndot = A.y(1,spot);
            pressure.cal = A.y(2,spot);
            bed.ini = A.x(1,spot);
            decompose.lamda = (decompose.ndot / ndot.initial.h2o2) * 100;

    case 1
            decompose.lamda = (decompose.ndot / ndot.initial.h2o2) * 100;
    
     % % mole fraction
       X.gas.o2 = ndot.o2 / (ndot.o2 +evap_cal * (ndot.h2o2 + ndot.h2o));
       X.gas.h2o2 = ndot.h2o2*evap_cal / (ndot.o2 + evap_cal*(ndot.h2o2 + ndot.h2o));
       X.gas.h2o = ndot.h2o*evap_cal / (ndot.o2 + evap_cal* (ndot.h2o2 + ndot.h2o));
       X.liquid.h2o2 = ndot.h2o2 / (ndot.h2o2 + ndot.h2o);
       X.liquid.h2o = ndot.h2o / (ndot.h2o2 + ndot.h2o);

    % % flow properties calculation

       ro = density(temp);
       volflux.liquid = (1-evap_cal) * ((ndot.h2o2 * mw.h2o2 / ro.h2o2) + (ndot.h2o * mw.h2o / ro.h2o));
       volflux.gas = 0.001 * (ndot.o2 + evap_cal * (ndot.h2o2 + ndot.h2o)) * (R * temp.cal.K / pressure.cal);
       ro.total = (ndot.h2o2 * mw.h2o2 + ndot.h2o * mw.h2o + ndot.o2 * mw.o2) / (volflux.liquid + volflux.gas);
       void = volflux.gas / (volflux.liquid + volflux.gas);

            % % Molar concentration of all species

            concen.h2o2 = ndot.h2o2 / (volflux.liquid + volflux.gas);
            concen.h2o = ndot.h2o / (volflux.liquid + volflux.gas);
            concen.o2 = ndot.o2 / (volflux.liquid + volflux.gas);

            % % Dynamic and kinematic viscosities of the homogeneous flow

            vis_inputs.R = R;
            vis_inputs.temp = temp;
            vis_inputs.ndot = ndot;
            vis_inputs.X = X;
            vis_inputs.mw = mw;
            
            mu = viscosity(vis_inputs);
            mu.mix = mu.liquid.mixture * (1-void) + (mu.gas.mixture * void);  % Dynamic viscosities
            dyv =  mu.mix / ro.total;                                         % kinematic viscosities
            Rep = supV * bed.pellet_dia / dyv;                                % Reynolds number

            % Assume that thee concentration of H2O2 on the catalyst surface is identical to that of the H2O2 in the flow %

            concen.surface.h2o2 = concen.h2o2;

            supV = (volflux.liquid + volflux.gas) / bed.area_eff;
            rdot = Af2 * (exp(-Ef2 / R * temp.cal.K)) * Ns *(K1 * concen.surface.h2o2 / (1 + K1 * concen.surface.h2o2));

            ode_inputs.volflux = volflux;
            ode_inputs.Rep = Rep;
            ode_inputs.supV = supV;
            ode_inputs.rdot = rdot;
            ode_inputs.bed = bed;
            ode_inputs.ro = ro;
            xspan = [bed.ini bed.ini + bed.dx/2 bed.ini + bed.dx];
            A = ode45(@(x,Y) fun_lamda2(x,Y,ode_inputs),xspan,[pressure.cal]);
            spot = length(A.x);

            pressure.cal = A.y(1,spot);
            bed.ini = A.x(1,spot);
            decompose.lamda = 1;

    end


    pressure_plot(i) = pressure.cal;
    decompose_ndot_plot(i) = decompose.ndot;
    decompose_lamda_plot(i) = decompose.lamda;
    temp_plot(i) = temp.cal.K;
    evap_plot(i) = evap_cal;

    printsol.temp = temp;
    printsol.evap = evap_cal;
    printsol.decompose.lamda = decompose.lamda;

    display = [printsol.evap, printsol.temp.cal.K, printsol.decompose.lamda];
    fprintf('evaporation = %6.4f where temperature = %8.5f\n and decompose rate = %8.5f \n',display)

end

    
    
    

       
       















