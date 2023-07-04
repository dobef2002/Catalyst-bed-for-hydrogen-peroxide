function boil_temp = boiling(WT_h202 , pressure)

wt = WT_h2o2;     %weight concentration of h2o2 in %

P = pressure.cal * 0.000001; % pa to Mpa

boil_temp.K = (412.65 +101.03 * wt) + (91.163 + 21.567 * wt) * log10(P) + (23.348 + 5.0267 * wt)*(log10(P))^2 + (5.4631 + 0.9383 * wt) * (log10(P))^3;

boil_temp.C = boil_temp.K - 273.15;

end