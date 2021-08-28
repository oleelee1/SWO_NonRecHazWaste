function MolarRatio = MROrg(C,H,N,O)
% function to determine the molar ratio of the oxygen required for
% oxidation
A = [0 -1 0 0;
     0 0 -2 0;
     0 0 0 -2;
     2 -2 -1 0];
b = [-C; -H; -N; -O];
X = A\b;
MolarRatio = X(1);