function feedpumppower = PumpPow(P1,P2,eff)
%function to claculate the pumping power required for feed pump
WaterDataNIST = readtable('Water_25degC_1to250bar');
P = WaterDataNIST(:,2);
Dens = WaterDataNIST(:,3);
Enth = WaterDataNIST(:,6);
P = table2array(P);
Dens = table2array(Dens);
Enth = table2array(Enth);

[P1diff, P1ind] = min(abs(P-P1));
[P2diff, P2ind] = min(abs(P-P2));

sumdens = sum(Dens(P1ind:P2ind));
AvDens = sumdens/length(P1ind:P2ind);

work = (100000*(P2-P1))/AvDens; %J/kg
feedpumppower = work/eff;
end

