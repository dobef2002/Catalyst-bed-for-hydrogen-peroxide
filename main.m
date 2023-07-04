

%%% catalytic bed calculation  %%% ode45
%%% unit in m-kg-s
%%% Calculate the decomposition rate, pressure drop and temperature in a catalystic bed for hydrogen peroxide.
%%% Use a simplfied model to demonstrate the realationship between evapration, enthalpy conservation and flow pressure for hydrogen peroxide decomposition activity in a catalyst bed.
%%% The program written with reference to the following article.
%%% (1) A.pasini et al"Reduced-Order Model for H2O2 Catalytic Reactor Performance Analysis"
%%% (2) "Scaling of catalyst bed for hydrogen peroxide monopropellant thrusters using catalytic decomposition modeling"

%%% solve decompose rate and pressure drop ODE with ODE45

clc;
clear;
global volflux Rep supV rdot bed.porosity ro bed checkflag.div

%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
user_input.WT = 0.90;  % weight percentage
user_input.bed.dia = 0.029; %(m)
user_input.bed.length = 0.0294; %(m)
user_input.pellet_dia = 0.0018;  %(m)
user_input.bed.segment = 10000;  %  bed length divide
user_input.mdot.total = 0.067; %kg/s
user_input.temp.enviroment.K = 300 ; % (k)
user_input.pressure.inlet = 2550000; %(pa)  initial pressure
user_input.Ns = 1;
user_input.K1 = 1;
user_input.Aspecific = 19.7941;
user_input.Ef2 = 15000; %kj/kmol
%%%%%%%%%%%%%%%%%%%%%%%%% USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % propellant properties
R = 8.314; % kJ/(kmol*k) ideal gas constant
mw.h2o2 = 34;
mw.h2o = 18;
mw.o2 = 32;
WT = user_input.WT; % weight percentage



% % bed dimension
bed.dia = user_input.bed.dia;             % catalyst bed diameter (m)
bed.length = user_input.bed.length ;        % catalyst bed length   (m)
bed.pellet_dia = user_input.pellet_dia;     % catalyst pellet or mesh diameter/wire diameter (m)
bed.porosity = 0.375+0.34*bed.pellet_dia/bed.dia;    % bey equation
bed.segment = user_input.bed.segment;
bed.dx = bed.length/bed.segment;
bed.area = (pi*bed.dia^2)/;  %(m^2)
bed.ini = 0;

% % catalyst and reaction rate properties
reac.Ns= user_input.Ns; reac.K1 = user_input.K1;
reac.Aspecific = user_input.Aspecific;                   % Arrhenius pre-exponential factor per unit specific area,[kmol/m2⋅s]
reac.Ef2 =  user_input.Ef2;                           % activation energy, kJ/kmol
reac.as = (1-bed.porosity)6/bed.pellet_dia;  % specific area of catalyst, m^(-1)
reac.Af2 = reac.as*reac.Aspecific;                     % Arrhenius pre-exponential factor, [kmol/m3⋅s]
bed.area_eff = bed.area*(1-bed.porosity)     % effective area

% % mass flow rate kg/s
mdot.total = user_input.mdot.total;                         %
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

temp.enviroment.K = user_input.temp.enviroment.K; %(k)
temp.enviroment.C = temp.enviroment.K - 273.15; % (C)
temp.dTp = 0.001;
temp.initial = temp.enviroment.K + temp.dTp;
temp.cal.K = temp.initial;
temp.cal.C = temp.cal.K - 273.15
temp.ref = 298.15;
h = enthalpy(temp);
h.total.initial = ndot.initial.h2o2*h.enthalpy.liquid.h2o2 + ndot.initial.h2o*h.enthalpy.liquid.h2o;

% % pressure
pessure.inlet = user_input.pressure.inlet; % pa
pressure.drop = 0;
pressure.cal = pessure.inlet - pressure.drop;


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
decompose_ndot_plot = zeros(bed.segment,1);
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
elseif (ndot.initial.h2o2 - decompose.ndot) < 0.00000000000000000000000000000000000000000001
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
       volflux.gas = 1000 * (ndot.o2 + evap_cal * (ndot.h2o2 + ndot.h2o)) * (R * temp.cal.K / pressure.cal);
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


            % Assume that thee concentration of H2O2 on the catalyst surface is identical to that of the H2O2 in the flow %

            concen.surface.h2o2 = concen.h2o2;

            supV = (volflux.liquid + volflux.gas) / bed.area_eff;
            rdot = Af2 * (exp(-reac.Ef2 / R * temp.cal.K)) * reac.Ns *(reac.K1 * concen.surface.h2o2 / (1 + reac.K1 * concen.surface.h2o2));
            Rep = supV * bed.pellet_dia / dyv;                                % Reynolds number

            ode_inputs.volflux = volflux;
            ode_inputs.Rep = Rep;
            ode_inputs.supV = supV;
            ode_inputs.rdot = rdot;
            ode_inputs.bed = bed;
            ode_inputs.ro = ro;

            xspan = [bed.ini bed.ini + bed.dx/2 bed.ini + bed.dx];
            A = ode45(@(x,y) fun_dlamda(x,y,ode_inputs),xspan,[decompose.ndot,pressure.cal]);
            %A = ode45(@fun_dlamda,xspan,[decompose.ndot,pressure.cal]);
            spot = length(A.x);
            decompose.ndot = A.y(1,spot);
            pressure.cal = A.y(2,spot);
            bed.ini = A.x(1,spot);
            decompose.lamda = (decompose.ndot / ndot.initial.h2o2) * 100;

    case 1
            decompose.lamda = (decompose.ndot / ndot.initial.h2o2) * 100;
            evap_cal = 1;

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

     % % mole fraction
       X.gas.o2 = ndot.o2 / (ndot.o2 +evap_cal * (ndot.h2o2 + ndot.h2o));
       X.gas.h2o2 = ndot.h2o2*evap_cal / (ndot.o2 + evap_cal*(ndot.h2o2 + ndot.h2o));
       X.gas.h2o = ndot.h2o*evap_cal / (ndot.o2 + evap_cal* (ndot.h2o2 + ndot.h2o));
       X.liquid.h2o2 = ndot.h2o2 / (ndot.h2o2 + ndot.h2o);
       X.liquid.h2o = ndot.h2o / (ndot.h2o2 + ndot.h2o);

    % % flow properties calculation

       ro = density(temp);
       volflux.liquid = 0;
       % volflux.liquid = (1-evap_cal) * ((ndot.h2o2 * mw.h2o2 / ro.h2o2) + (ndot.h2o * mw.h2o / ro.h2o));
       volflux.gas = 1000 * (ndot.o2 + evap_cal * (ndot.h2o2 + ndot.h2o)) * (R * temp.cal.K / pressure.cal);
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


            % Assume that thee concentration of H2O2 on the catalyst surface is identical to that of the H2O2 in the flow %

            concen.surface.h2o2 = concen.h2o2;

            supV = (volflux.liquid + volflux.gas) / bed.area_eff;
            rdot = reac.Af2 * (exp(-reac.Ef2 / R * temp.cal.K)) * reac.Ns *(reac.K1 * concen.surface.h2o2 / (1 + reac.K1 * concen.surface.h2o2));
            Rep = supV * bed.pellet_dia / dyv;                                % Reynolds number

            ode_inputs.volflux = volflux;
            ode_inputs.Rep = Rep;
            ode_inputs.supV = supV;
            ode_inputs.rdot = rdot;
            ode_inputs.bed = bed;
            ode_inputs.ro = ro;
            xspan = [bed.ini bed.ini + bed.dx/2 bed.ini + bed.dx];
            A = ode45(@(x,y) fun_lamda2(x,y,ode_inputs),xspan,[pressure.cal]);
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