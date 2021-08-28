clear;
clc;
close all
format longg
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')
warning('OFF', 'MATLAB:legend:IgnoringExtraEntries')
    
%% Parameters 
T0 = 25;
T1 = [400]%, 450, 475, 500];
T20 = 800;
T3 = T1+10;
TS1 = 25;
TS2 = 400;

X = 0.597

RTwant = 20;

mass_basis = 1000; %kg/hr
SR = 1;
wt = 0.002415;

%% Create arrays for each variable
lengthP1 = length(T1);

RT_P1 = ones(lengthP1,1);

RT_P1_X = ones(lengthP1,length(X));

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

                    Massflow_Mea = mass_basis*wt; %kg/hr
                    Molflow_Wat = mass_basis*(1-wt)/Mrwat; %kmol/hr
                    Molflow_Mea = Massflow_Mea/MrMea; %kmol/hr
                    Molflow_Ox = Molflow_Mea*SR; % kmol/hr

                    Massflow_Ox = Molflow_Ox*Mrox; %kg/hr
                    Massflow_wat = mass_basis - Massflow_Mea - Massflow_Ox;

                    % Calculate enthalpy required 

                    TotalEnergy = H_inlet_req*mass_basis; %kJ/hr  - (+ H_reac)

                    % Calculate the Moles Required to be Oxidised To Provide This heat 

                    C = 2; H = 7; O = 1;
                    Delhr = DelHR_mass(C,H,O,MrMea); % Kj/kgMea
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

             %% Attempt at integrating the 1/k expression to calculate residence time

                    a = 4533.6; %MEA
                    b = 3.6297; %MEA
                    %a = 8558; %MEA IPA
                    %b = 8.8797; %MEA IPA
                    %a = 6294.5; %MEA PG
                    %b = 6.3787; %MEA PG
                    %a = 11462; %3MP
                    %b = 11.7; %3MP
                    %a = 8290.1; %3MP IPA
                    %b = 7.4145; %3MP IPA
                    %a = 8275.5; %3MP PG
                    %b = 8.026; %3MP PG
                    
                    Xa = 0:0.001:X(x);
                    c = 1-Xa;
                    d = (T1(P1)+273) + ((wt.*Xa.*Delhr)./(CpAv));
                    e = 1./d;
                    f = ((-a).*e)+b;
                    g = exp(f);
                    h = c.*g;
                    F = 1./h;
                    S = trapz(Xa,F);
                    RT = double(S);
                    
                    T2(i) = (T1(P1)+273) + ((wt*X(x)*Delhr)/(CpAv));
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

RT_P1(P1) = RT;
    end
RT_P1_X(:,x) = RT_P1
[diffval indRT] = min(abs(RT_P1_X-RTwant))
[X(indRT), RT_P1_X(indRT)]
end
