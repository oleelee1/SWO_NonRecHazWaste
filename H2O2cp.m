function H2O2_specific_heat = H2O2cp(T)
% function to calculate the specific heats of hydrogen peroxide 

A	=34.25667;
B	=55.18445;
C	=(-35.15443);
D	=9.087440;
E	=(-0.422157);
F	=(-149.9098);
G	=257.0604;
H	=(-136.1064);
t = T./1000; %kelvin

cp = A + B.*t + C.*(t.^2) + D.*(t.^3) + E./(t.^2);
H2O2_specific_heat = cp; %j/mol.K
end
