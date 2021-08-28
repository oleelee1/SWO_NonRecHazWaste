function SteamPower_T1 = ExtrapSteam(TS2,Pi,Po,mS,T1,X)
% function to work out the power generated from the steam turbine given

mS = mS./3600; % kg/s


%% Steam data - T - Temp; E - Entropy; h - Enthalpy

if Pi == 30
    steamdata_inlet = readtable('S_30');
elseif Pi == 45
        steamdata_inlet = readtable('S_45');
elseif Pi == 17.5 
    steamdata_inlet = readtable('S_17_5');
end

if Po == 0.07
    steamdata_outlet = readtable('S_0_07');
elseif Po == 0.3
    steamdata_outlet = readtable('S_03');
elseif Po == 1
    steamdata_outlet = readtable('S_1');
elseif Po == 3
    steamdata_outlet = readtable('S_3');
elseif Po == 10
    steamdata_outlet = readtable('S_10');
end

T_inlet = steamdata_inlet(:,1);
T_outlet = steamdata_outlet(:,1);
T_inlet = table2array(T_inlet);
T_outlet = table2array(T_outlet);

E_inlet = steamdata_inlet(:,7);
E_outlet = steamdata_outlet(:,7);
E_inlet = table2array(E_inlet);
E_outlet = table2array(E_outlet);

h_inlet = steamdata_inlet(:,6);
h_outlet = steamdata_outlet(:,6);
h_inlet = table2array(h_inlet);
h_outlet = table2array(h_outlet);

[Tin, Tindin] = min(abs(TS2-T_inlet));

Entropy = E_inlet(Tindin); % inlet entholpy = outlet entropy (assume isentropic)

[Eout, Eindout] = min(abs(Entropy - E_outlet));

ho = h_outlet(Eindout);
hi = h_inlet(Tindin);
Q = 0;


%% plotting power 

W = Q-mS*(ho-hi);
W = W*3600; % KJ/hr
Tout = T_outlet(Eindout)

hp = W/2684.5195;

%% extrapolation
horse = [5, 30, 175, 325, 425];
effic = [0.15, 0.2, 0.25, 0.3, 0.4];

coefficients = polyfit(horse, effic, 3);
EffExtrap = polyval(coefficients, hp);
   
if EffExtrap >= 0.4
    EffExtrap = 0.4;
end

SteamPower_T1 = W*EffExtrap/3600; %KW

    disp('horsepower')
    disp(hp)
    %disp('efficiency')
    %%disp(eff)
    %disp('Power Produced (KW)')
    %disp(SteamPower_T1)
    %disp('Temperature of the outlet steam (deg C)')
    %disp(Tout)
end