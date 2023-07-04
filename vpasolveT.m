function outputs = vpasolveT(inputs)

          temp = inputs.temp;
          ndot = inputs.ndot;
          pressure = inputs.pressure;
          h = inputs.h;
          checkflag = inputs.checkflag;

switch checkflag.evap
       case 0
        syms evap T

             X.gas.o2 = ndot.o2 / (ndot.o2 +evap * (ndot.h2o2 + ndot.h2o));
          %  X.gas.h2o2 = ndot.h2o2*evap_cal / (ndot.o2 + evap_cal*(ndot.h2o2 + ndot.h2o));
          %  X.gas.h2o = ndot.h2o*evap_cal / (ndot.o2 + evap_cal* (ndot.h2o2 + ndot.h2o));
             X.liquid.h2o2 = ndot.h2o2 / (ndot.h2o2 + ndot.h2o);
             X.liquid.h2o = ndot.h2o / (ndot.h2o2 + ndot.h2o);

             % % partial  pressure

             pressure = p_saturationT(temp.pressure,T);
             pressure.partial.h2o2 = pressure.sat.h2o2 * X.liquid.h2o2;
             pressure.partial.h2o = pressure.sat.h2o * X.liquid.h2o;
             pressure.partial.o2 = pressure.cal * X.gas.o2;

             % % determine temperature and evaporation rate

             eqn1 = pressure.cal - (pressure.partial.h2o2 + pressure.partial.h2o + pressure.partial.o2) == 0;

                 % % enthalpy conservation

                 h = enthalpyT(temp,h,T);
                 h.gas = (ndot.o2 * h.enthalpy.gas.o2 + evap * (ndot.h2o2 * h.enthalpy.gas.h2o2 + ndot.h2o * h.enthalpy.gas.h2o));
                 h.liquid = (1-evap) * (ndot.h2o2 * h.enthalpy.liquid.h2o2 + ndot.h2o * h.enthalpy.liquid.h2o);

             eqn2 = h.total.initial - h.gas - h.liquid == 0;


             assume(T,'real')
             assume(evap, 'real')
             T_cal = temp.cal.K;
             solx = vpasolve([eqn1,eqn2],[evap,T],[evap_cal,T_cal]); % use vpasolve to solve unknown evap and T
             evap_sol = round(double(solx.evap),20);
             T_sol = round(double(solx.T),20);

             [evap_cloum,evap_row] = size(evap_sol);
             [T_cloum,T_row] = size(T_sol);

             for ii = 1:T_cloum
                    if evap_sol(ii) >= 0 && evap_sol(ii) < 1
                          if T_sol(ii) < 299
                             res(ii) = 10000;
                             continue
                          elseif T_sol(ii) > 1600
                             res(ii) = 10000;
                              continue
                          else
                              res(ii) = abs(real(T_sol(ii))-T_cal);
                          end
                    elseif  evap_sol(ii) >= 1
                          if T_sol(ii) < 299
                             res(ii) = 10000;
                             continue
                          elseif T_sol(ii) > 1600
                              res(ii) = 10000;
                              continue
                          else
                              res(ii) = abs(real(T_sol(ii))-T_cal);
                          end

                          continue
                     end
              end

           [res_number,res_locate] = max(res);
           T_cal = T_sol(res_locate);
           evap_cal = evap_sol(res_locate);
           temp.cal.K = T_cal;

        case 1
            evap_cal = 1 ;
            evap = 1;
            syms T2

            % % enthalpy conservation

            h = enthalpyT(temp,h,T2);
            h.gas = (ndot.o2 * h.enthalpy.gas.o2 + evap * (ndot.h2o2 * h.enthalpy.gas.h2o2 + ndot.h2o * h.enthalpy.gas.h2o));
            h.liquid = (1-evap)*(ndot.h2o2 * h.entalpy.liquid.h2o2 + ndot.h2o * h.entalpy.liquid.h2o);
            eqn2 = h.total.initial - h.gas - h.liquid == 0;
            T_cal = temp.cal.K;
            assume(T2,'real')
            solx = vpasolve (eqn2,T2,T_cal);
            T_sol = round(double(solx),20);
            [T_cloum,T_row] = size(T_sol);

            count = 0;

            for ii = 1:T_cloum

                          if T_sol(ii) < 299
                              res(ii) = 10000;
                             continue
                          elseif T_sol(ii) > 1600
                              res(ii) = 10000;
                              continue
                          else
                              res(ii) = abs(real(T_sol(ii))-T_cal);
                              count = count +1;
                          end

            end

            if count == 0
            checkflag.div = 1;
            disp('program has no solution')
            end

            %[res_number,res_locate] = max(res);
            [res_number,res_locate] = min(res);
            T_cal = T_sol(res_locate);
            evap_cal = 1;
            temp.cal.K = T_cal;
            temp.cal.K = real(T_cal(1));
        end
end



