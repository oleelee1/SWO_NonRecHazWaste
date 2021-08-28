clear;
clc;
close all
format longg
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')
warning('OFF', 'MATLAB:legend:IgnoringExtraEntries')
    
%% Parameters 
T0 = 25;
%T1 = 350:10:500;
T1 = linspace(350,500,20);
T20 = 800;
T3 = T1+10;
TS1 = 25;
TS2 = 400;
X = [0.8, 0.85, 0.9, 0.95 ,0.99];
mass_basis = 1000; %kg/hr
SR = 1;
wt = 0.05;

%% Create arrays for each variable

T2_P1 = ones(length(T1),1);
T1_P1 = ones(length(T1),1);
CpAv_P1 = ones(length(T1),1);
RT_P1 = ones(length(T1),1);
C0f_P1 = ones(length(T1),1);
Massflow_Mea_P1 = ones(length(T1),1);
mS_P1 = ones(length(T1),1);
ReactorVol_P1 = ones(length(T1),1);
ReactorLength_P1 = ones(length(T1),1);
power_water_P1 = ones(length(T1),1);
power_organic_P1 = ones(length(T1),1);
hex_steam_A_L_P1 = ones(length(T1),3);
hex_feed_A_L_P1 = ones(length(T1),3);
wtf_P1 = ones(length(T1),1);
H_steam_avail_P1 = ones(length(T1),1);
CompWork_HP_P1 = ones(length(T1),1);
PumpWork_HP_P1 = ones(length(T1),1);

T2_P1_X = ones(length(T1),length(X));
RT_P1_X = ones(length(T1),length(X));
T1_P1_X = ones(length(T1),length(X));
CpAv_P1_X = ones(length(T1),length(X));
C0f_P1_X = ones(length(T1),length(X));
Massflow_Mea_P1_X = ones(length(T1),length(X));
mS_P1_X = ones(length(T1),length(X));
ReactorVol_P1_X = ones(length(T1),length(X));
ReactorLength_P1_X = ones(length(T1),length(X));
power_water_P1_X = ones(length(T1),length(X));
power_organic_P1_X = ones(length(T1),length(X));
hex_steam_A_L_P1_X = ones(length(T1),length(X)*3);
hex_feed_A_L_P1_X = ones(length(T1),length(X)*3);
wtf_P1_X = ones(length(T1),length(X));
H_steam_avail_P1_X = ones(length(T1),length(X));
CompWork_HP_P1_X = ones(length(T1),length(X));
PumpWork_HP_P1_X = ones(length(T1),length(X));

%% load tables 
WaterDataNIST = readtable('ResearchProject_AccCpData_kg');


%% Start for loop to change paramaters (T1, mS, X, etc)

for x = 1:length(X)
    for P1 = 1:length(T1)
        
        Vmat = Vol_m(mass_basis,T1(P1));
        V = Vmat(1);
        Dens = Vmat(2);

            %% Start iteration and Calculate enthalpy of each stream

            C0 = ones(2,1);
            T2 = ones(2,1);

            i = 2;
            T2(i) = T20;
            T2(i-1) = 0;
            
            while abs(T2(i)-T2(i-1))>0.5
                
                T2(i-1) = T2(i);
                
                    % NIST data
                    T = WaterDataNIST(:,1); % deg C
                    Enth = WaterDataNIST(:,6); % kj/kg

                    T = table2array(T);
                    Enth = table2array(Enth);
                    
                    TInt = T(470:601);
                    EnthInt = Enth(470:601);

                    TExtrap = [1000:2000]';

                    EnthExtrap = interp1(TInt,EnthInt,1000:2000,'linear','extrap');
                    EnthExtrap = EnthExtrap';

                    T = [T;TExtrap];
                    Enth = [Enth;EnthExtrap];

                    [Tval0, ind0] = min(abs(T-T0));
                    [Tval1, ind1] = min(abs(T-T1(P1)));
                    [Tval2, ind2] = min(abs(T-T2(i-1)));
                    [Tval3, ind3] = min(abs(T-T3(P1)));
                    [TvalS1, indS1] = min(abs(T-TS1));
                    [TvalS2, indS2] = min(abs(T-TS2));

                    Enth0 = Enth(ind0);
                    Enth1 = Enth(ind1);
                    Enth2 = Enth(ind2);
                    Enth3 = Enth(ind3);
                    EnthS1 = Enth(indS1);
                    EnthS2 = Enth(indS2); %kj/kg

                    H_inlet_req = Enth1 - Enth0; % KJ/kg 
                    Enth4 = Enth3 - Enth1 + Enth0;
                    
                    [Enth4val, ind4] = min(abs(Enth-Enth4));
                    T4 = T(ind4);

                    H_inlet_avail = Enth3 - Enth4; % KJ/kg
                    H_reac = Enth2 - Enth1;

                    H_steam_avail = Enth2 - Enth3; %Kj/kg
                    H_steam_req = EnthS2 - EnthS1;
                    
            %% start wt iteration

                    % Calculate mass

                    C = 2; H = 7; N = 1; O = 1;
                    MROx = MROrg(C,H,N,O); 

                    Mrwat = 18;
                    MrMea = 61;
                    Mrox = 34;
                    MrIPA = 60.1;
                    MrPG = 76.09;


                    Massflow_Mea = mass_basis*wt; %kg/hr
                    Molflow_Wat = mass_basis*(1-wt)/Mrwat; %kmol/hr
                    Molflow_Mea = Massflow_Mea/MrMea; %kmol/hr
                    Molflow_Ox = Molflow_Mea*SR; % kmol/hr
                    Molflow_IPA = Molflow_Ox;
                    Molflow_PG = Molflow_Ox;

                    Massflow_Ox = Molflow_Ox*Mrox; %kg/hr
                    Massflow_IPA = Molflow_IPA*MrIPA; %kg/hr
                    Massflow_PG = Molflow_PG*MrPG; %kg/hr
                                        
                    Massflow_wat = mass_basis - Massflow_Mea - Massflow_Ox - Massflow_IPA;

                    % Calculate enthalpy required 

                    TotalEnergy = H_inlet_req*mass_basis; %kJ/hr  - (+ H_reac)

                    % Calculate the Moles Required to be Oxidised To Provide This heat 

                    C = 2; H = 7; O = 1;
                    Delhr = DelHR_mass(C,H,O,MrMea); % Kj/kgMea
                    DelhrIPA = DelHR_mass(3,8,1,MrIPA); %KJ/kg
                    DelhrPG = DelHR_mass(3,8,2,MrPG);
                    
                    %kgMea = TotalEnergy/(Delhr*X(x)); %kgMea/hr
                    %wt(i) = kgMea/mass_basis(mb); %kgMea/kgfeed
                    kgMea = Massflow_Mea;

                    kgCOD = (Massflow_Mea/MrMea)*MROx*32; % kgOxygen/hr
                    wtCOD = kgCOD/mass_basis; % kgCOD/kgFeed

                    C0kmolkg = wt/MrMea; % kmol/kgfeed
                    C0kmolvol = C0kmolkg*mass_basis/V; %kmol/m3
                    C0molvol = C0kmolvol*1000; %mol/m3
                    CODmolvol = (C0molvol)*MROx*32;

                    C0gvol = C0molvol*MrMea; %g/m3
                    C0 = C0gvol/1000; %kg/m3

                    kgCODf = kgCOD;
                    wtCODf = wtCOD;
                    C0f = C0molvol;
                    
                    mS = (H_steam_avail*mass_basis)/H_steam_req;
                    
                    H_steam_rem = (H_steam_avail*mass_basis)-(H_steam_req*mS);

                    % Calculate Q             
                    Qf_avail = 1000*H_inlet_avail.*mass_basis./3600; % J/s
                    Qf_req = 1000*H_inlet_req.*mass_basis./3600; % J/s
                    Qs = 1000*H_steam_req.*mS./3600; %J/s
                    
           %% Calculate pump and compressor power requirements
           
           eff = 0.72;
           specwork = PumpPow(1,250,eff); %J/kg
           PumpWork_KW = specwork * (Massflow_Mea+Massflow_wat)/(3600*1000); %KJ/s or KW
           PumpWork_HP = PumpWork_KW*1.341;
           
           eff = 0.82;
           T1comp = T0+273;
           specwork = H2O2comp(T1comp,1,250,eff); %J/kg
           CompWork_KW = specwork * Massflow_Ox/(3600*1000); %KW
           CompWork_HP = CompWork_KW * 1.341;
           

            %% Calculate average Cp in reactor;
  
                    T = WaterDataNIST(:,1); % deg C
                    Cp = WaterDataNIST(:,9); % kj/kg.K

                    T = table2array(T);
                    Cp = table2array(Cp);   
                    
                    TInt = T(470:601);
                    CpInt = Cp(470:601);
                    
                    TExtrap = [1000:2000]';

                    coefficients = polyfit(TInt, CpInt, 2);
                    CpExtrap = polyval(coefficients, 1000:2000);
                    CpExtrap = CpExtrap';

                    T = [T;TExtrap];
                    Cp = [Cp;CpExtrap];

                    [T1diff, T1ind] = min(abs(T-T1(P1)));
                    [T2diff, T2ind] = min(abs(T-T2(i-1)));
                   
                    CpSum = sum(Cp(T1ind:T2ind));
                    CpAv = CpSum/length(Cp(T1ind:T2ind)); %Kj/kg.k
                    %CpAv = 3.5;

            %% Attempt at integrating the 1/k expression to calculate residence time

                    %a = 4533.6; %MEA
                    %b = 3.6297; %MEA
                    a = 8558; %MEA IPA
                    b = 8.8797; %MEA IPA
                    %a = 6294.5; %MEA PG
                    %b = 6.3787; %MEA PG
                    %a = 11462; %3MP
                    %b = 11.7; %3MP
                    %a = 8290.1; %3MP IPA
                    %b = 7.4145; %3MP IPA
                    %a = 8275.5; %3MP PG
                    %b = 8.026; %3MP PG
                    
                    Xa = 0:0.01:X(x);
                    c = 1-Xa;
                    d = (T1(P1)+273) + ((wt.*Xa.*Delhr)./(CpAv)) + ((Xa.*Massflow_IPA*DelhrIPA)/(mass_basis*CpAv));
                    e = 1./d;
                    f = ((-a).*e)+b;
                    g = exp(f);
                    h = c.*g;
                    F = 1./h;
                    S = trapz(Xa,F);
                    RT = double(S);
                    
                    T2(i) = (T1(P1)+273) + ((wt*X(x)*Delhr)/(CpAv))+ ((X(x).*Massflow_IPA*DelhrIPA)/(mass_basis*CpAv));
                    T2(i) = T2(i)-273;

            end
             T2f = T2(i(end));
                        
            
%% Calculate the Volume and length of reactor required

VFlow = V/(1*60*60);
% RT = V/(VFlow*A)
% V = A*L
% A = Pi*R^2
% RT = L/VFlow

R = 0.0508;
SA = 3.14159*R^2;

ReactorVol = RT*VFlow;
ReactorLength = ReactorVol/SA;

ActReactorLength = 10;
ActReactorVol = ActReactorLength*SA;
ActRT = ActReactorVol/VFlow;

%% Reynolds 

u=VFlow/((pi*(R^2)));
mu=30*(10^-6);
Re=(Dens*u*2*R)/mu;

%% Calculate the power produced by the steam turbine

Power = ExtrapSteam(TS2,30,0.07,mS,T1(P1),X(x));

%% Calculate the power produced per kg wastewater/organic treated

power_water = (Power*3600)/mass_basis;
power_organic = (Power*3600)/Massflow_Mea;

%% HEX

% Q = UAlmdT

% Feed heat exchanger
U = 3000;
D = 25e-3;
L = 7.32;
%L_Dshell = 10;
%D_shell = L/L_Dshell;


delT1 = T3(P1)-T1(P1);
delT2 = T4-T0;
lmdT = (delT1-delT2)/(log(delT1/delT2));
hex_feed_A_L = HEX(lmdT,Qf_req,U,D,L);
%hex_feed_A_L = HEX(lmdT,Q(end),U,D,L);


% Steam heat exchanger
U = 3000;
D = 16e-3;
L = 1.83;
L_Dshell = 10;
D_shell = L/L_Dshell;

delT1 = T2f-TS2;
delT2 = T3(P1)-TS1;
lmdT = (delT1-delT2)/(log(delT1/delT2));
hex_steam_A_L = HEX(lmdT,Qs,U,D,L);
    
%% Results
           
            T2_P1(P1) = T2f;
            T1_P1(P1) = T1(P1);
            CpAv_P1(P1) = CpAv;
            Massflow_Mea_P1(P1) = Massflow_Mea;
            C0f_P1(P1) = C0f;
            mS_P1(P1) = mS;
            RT_P1(P1) = RT;
            ReactorVol_P1(P1) = ReactorVol;
            ReactorLength_P1(P1) = ReactorLength;
            power_water_P1(P1) = power_water;
            power_organic_P1(P1) = power_organic;
            hex_steam_A_L_P1(P1,1:3) = hex_steam_A_L;
            hex_feed_A_L_P1(P1,1:3) = hex_feed_A_L;
            wtf_P1(P1) = wt;
            H_steam_avail_P1(P1) = H_steam_avail;
            CompWork_HP_P1(P1) = CompWork_HP;
            PumpWork_HP_P1(P1) = PumpWork_HP;
            
            
    end    
    T2_P1_X(:,x) = T2_P1;
    T1_P1_X(:,x) = T1_P1;
    CpAv_P1_X(:,x) = CpAv_P1;
    RT_P1_X(:,x) = RT_P1;
    C0f_P1_X(:,x) = C0f_P1;
    Massflow_Mea_P1_X(:,x) = Massflow_Mea_P1; % rows - T2, columns - X
    mS_P1_X(:,x) = mS_P1;
    ReactorVol_P1_X(:,x) = ReactorVol_P1;
    ReactorLength_P1_X(:,x) = ReactorLength_P1;
    power_water_P1_X(:,x) = power_water_P1;
    power_organic_P1_X(:,x) = power_organic_P1;
    hex_steam_A_L_P1_X(:,(x*3-2):(x*3)) = hex_steam_A_L_P1;
    hex_feed_A_L_P1_X(:,(x*3-2):(x*3)) = hex_feed_A_L_P1;
    wtf_P1_X(:,x) = wtf_P1;
    H_steam_avail_P1_X(:,x) = H_steam_avail_P1;
    CompWork_HP_P1_X(:,x) = CompWork_HP_P1;
    PumpWork_HP_P1_X(:,x) = PumpWork_HP_P1;
    
end
    
        
        
        
%% Arrays

%disp('Concentration Required at different T2s and Conversions')
%disp(C0f_P1_X)
%disp('Mass fraction Required at different T2s and Conversions')
%disp(wtf_P1_X)
disp('Mass of Steam at different mass basis and Conversions')
disp(mS_P1)
%disp('Reactor Outlet Temperature')
%disp(T2_P1_X)
disp('Residence Time at different weight fractions and Conversions')
disp(RT_P1)
%disp('Reactor Volume at different T2s and Conversions')
%disp(ReactorVol_P1_X)
disp('Reactor Length at different weight fractions and Conversions')
disp(ReactorLength_P1)
%disp('Power produced per unit mass wastewater treated (kJ/kg)')
%disp(power_water_P1_X)
disp('Power produced per unit mass organic treated (kJ/kg)')
disp(power_organic_P1)
disp('Reactor Outlet Temperatures (deg C)')
disp(T2_P1)
%disp('Steam heat exchanger area, area per tube and nuwter of tubes required')
%disp(hex_steam_A_L_P1_X)
%disp('Feed heat exchanger area, area per tube and nuwter of tubes required')
%disp(hex_feed_A_L_P1_X)


%% plots 

figure(1)
i = X;
n = T1_P1_X+273;
m = T2_P1_X+273;

for i = 1:length(i)
    plot(n(:,i),m(:,i))
    hold on
end
hold off
grid on
xlabel('Reactor Inlet Temperature (K)')
ylabel('Reactor Outlet Temperature (K)')
legend('80', '85', '90', '95', '99')
leg = legend('show');
title(leg,'Conversion (%)')

figure(2)
i = X;
n = T1_P1_X + 273;

m = power_water_P1_X;

for i = 1:length(i)
    plot(n(:,i),m(:,i))
    hold on
end
hold off
grid on
xlabel('Reactor Inlet Temperature (K)')
ylabel('Power Generated per unit Wastewater (kW/kg)')

legend('80', '85', '90', '95', '99')
leg = legend('show');
title(leg,'Conversion (%)')


figure(3)
i = X;
n = T1_P1_X+273;
m = power_organic_P1_X;

for i = 1:length(i)
    plot(n(:,i),m(:,i))
    hold on
end
hold off
grid on
xlabel('Reactor Inlet Temperature (K)')
ylabel('Power Generated per unit Organic (kW/kg)')
legend('80', '85', '90', '95', '99')
leg = legend('show');
title(leg,'Conversion (%)')

figure(4)
i = X;
n = T1_P1_X + 273;

m = RT_P1_X;

for i = 1:length(i)
    plot(n(:,i),m(:,i))
    hold on
end
hold off
grid on
xlabel('Reactor Inlet Temperature (K)')
ylabel('Theoretical Residence Time (s)')

legend('80', '85', '90', '95', '99')
leg = legend('show');
title(leg,'Conversion (%)')

savefig(1,'T1_vs_T2_MEA_IPA_steam')
savefig(2,'T1_vs_poweroutwat_MEA_IPA_steam')
savefig(3,'T1_vs_powerout_MEA_IPA_steam')
savefig(4,'T1_vs_restime_MEA_IPA_steam')
save('T1_MEA_IPA_steam')