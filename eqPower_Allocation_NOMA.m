function [p1, p2] = eqPower_Allocation_NOMA(G1, G2, P, N)

syms P1 P2 
eqn1 = P1 < P2;
eqn2 = P1 > 0.5;
eqn3 = P2 > 0.5;
eqn4 = [G1*P1/N == G2*P2/(G2*P1+N), P == P1 + P2];
eqns = [eqn1 eqn2 eqn3 eqn4];
S = solve(eqns,[P1 P2], 'Real',true);

p1 = S.P1;
p2 = S.P2;
