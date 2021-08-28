function H2O2_compressor = H2O2comp(T1,P1,P2,eff)
% function to calculate the work required of the H2O2 compressor

array = H2O2cp(T1); 
cpox = array(1); % J/mol.K
cpox = cpox/34; % kJ/kg.K
cpox = cpox*1000; %J/kg.K
R = 1000*8.3145/34; %J/kg.K
cv = cpox-R;
gamma = cpox/cv;
specwork = (gamma/(gamma-1))*R*T1*(((P2/P1)^((gamma-1)/gamma))-1);
specwork = specwork/eff;
H2O2_compressor = specwork; % J/kg
end