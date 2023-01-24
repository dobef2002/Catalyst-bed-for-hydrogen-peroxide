function dYdx = fun_dlamda(x,Y,ode_inputs)

           volflux = ode_inputs.volflux;
           Rep = ode_inputs.Rep;
           supV = ode_inputs.supV;
           rdot = ode_inputs.rdot;
           bed = ode_inputs.bed;
           ro = ode_inputs.ro;


dYdx = zeros(2,1);
dYdx(1) = ((volflux.liquid + volflux.gas) / supV) * rdot; % d_decompose rate

% ergun equation
dYdx(2) = -(150*(1-bed.porosity) / Rep + 1.75) * (((1-bed.porosity) / bed.porosity^3) * ((ro.total * supV^2) / bed.pellet_dia));

% Tallmadge equation
%dYdx(2) = -(150*(1-bed.porosity) / Rep + (4.2 / (Rep^(1/6)))) * (((1-bed.porosity) / bed.porosity^3) * ((ro.total * supV^2) / bed.pellet_dia));

end