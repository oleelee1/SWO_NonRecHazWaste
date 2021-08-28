function HeatofReaction = DelHR_mass(C,H,O,Mr)
%Function to find the heat of reaction for a hydrocarbon based on C, H and
%O content
HeatofReaction = (415*C + 107*H - 193*O) * 1000 / Mr; % Kj/kgmea
end