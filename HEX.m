function hex_A_L = HEX(lmdT,Q,U,D,L)
% function to calculate the area of heat exchanger required for heat
% transfer

for i = 1:length(Q)

    A(i) = Q(i)/(lmdT*U);
    A_tube(i) = 3.14159*D*L;
    N_tube(i) = A(i)/A_tube(i);
    
end

hex_A_L = [A', A_tube', N_tube'];
end