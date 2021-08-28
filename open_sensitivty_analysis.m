clear;
clc;


Pdiff = open('Sens_5mM.mat')
vals = struct2cell(Pdiff)
vals_mat = cell2mat(vals)
sens_5mM = vals_mat';


Pdiff = open('Sens_10mM.mat')
vals = struct2cell(Pdiff)
vals_mat = cell2mat(vals)
sens_10mM = vals_mat';

Pdiff = open('Sens_15mM.mat')
vals = struct2cell(Pdiff)
vals_mat = cell2mat(vals)
sens_15mM = vals_mat';

Pdiff = open('Sens_20mM.mat')
vals = struct2cell(Pdiff)
vals_mat = cell2mat(vals)
sens_20mM = vals_mat';

Pdiff = open('Sens_425K.mat')
vals = struct2cell(Pdiff)
vals_mat = cell2mat(vals)
sens_425K = vals_mat';

Pdiff = open('Sens_450K.mat')
vals = struct2cell(Pdiff)
vals_mat = cell2mat(vals)
sens_450K = vals_mat';

Pdiff = open('Sens_475K.mat')
vals = struct2cell(Pdiff)
vals_mat = cell2mat(vals)
sens_475K = vals_mat';

Pdiff = open('Sens_500K.mat')
vals = struct2cell(Pdiff)
vals_mat = cell2mat(vals)
sens_500K = vals_mat';

ExpRT = [10;12;14;18;20];
format shortg
%Sensitivity_Analysis[1:5,1:9] = [ExpRT, sens_5mM(:,2), sens_5mM(:,1),...
    %ExpRT, sens_10mM(:,2), sens_10mM(:,1),...
    %ExpRT, sens_15mM(:,2), sens_15mM(:,1),...
    %ExpRT, sens_20mM(:,2), sens_20mM(:,1),...
    %ExpRT, sens_425K(:,2), sens_425K(:,1),...
    %ExpRT, sens_450K(:,2), sens_450K(:,1),...
    %ExpRT, sens_475K(:,2), sens_475K(:,1),...
    %ExpRT, sens_500K(:,2), sens_500K(:,1)]

Sensitivity_Analysis_1 = [ExpRT, sens_5mM(:,2), sens_5mM(:,1),...
    ExpRT, sens_10mM(:,2), sens_10mM(:,1)]
    
Sensitivity_Analysis_2= [ExpRT, sens_15mM(:,2), sens_15mM(:,1)...
    ExpRT, sens_20mM(:,2), sens_20mM(:,1)]
    
Sensitivity_Analysis_3 = [ExpRT, sens_425K(:,2), sens_425K(:,1),...
    ExpRT, sens_450K(:,2), sens_450K(:,1)]

Sensitivity_Analysis_4= [ExpRT, sens_475K(:,2), sens_475K(:,1),...
    ExpRT, sens_500K(:,2), sens_500K(:,1)]

sensanalysis = [Sensitivity_Analysis_1; Sensitivity_Analysis_2; Sensitivity_Analysis_3; Sensitivity_Analysis_4]
    
xlswrite('sensanalysis2',sensanalysis)



