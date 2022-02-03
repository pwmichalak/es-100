function [R0, R1, R2, Cap1, Cap2, T] = extract_2rc_params(A,B,C,D)
% return the resistor and capacitor values in a 2RC equivalent circuit
% battery model. Also return the transformation matrix T that converts
% between the obtained state space and the desired 2RC state space

% check that all matrices are of the correct size
assert(all(size(A)==[2 2]))
assert(all(size(B)==[2 1]))
assert(all(size(C)==[1 2]))
assert(all(size(D)==[1 1]))

% perform parameter extraction
R0 = D;

tau_plus = (-trace(A) + sqrt(trace(A)^2 - 4 * det(A))) / (2 * det(A));
tau_minus = (-trace(A) - sqrt(trace(A)^2 - 4 * det(A))) / (2 * det(A));
tau1 = max(tau_plus, tau_minus);
tau2 = min(tau_plus, tau_minus);

T11 = 1 / C(1) - (C(2) * A(1,1) * tau1 + C(2)) / ((C(1) * C(2) * A(1,1) - C(1)^2 * A(1,2)) * tau1 + C(1) * C(2));
T12 = 1 / C(1) - (C(2) * A(1,1) * tau2 + C(2)) / ((C(1) * C(2) * A(1,1) - C(1)^2 * A(1,2)) * tau2 + C(1) * C(2));
T21 = (A(1,1) * tau1 + 1) / ((C(2) * A(1,1) - C(1) * A(1,2)) * tau1 + C(2));
T22 = (A(1,1) * tau2 + 1) / ((C(2) * A(1,1) - C(1) * A(1,2)) * tau2 + C(2));
T = [T11 T12; T21 T22];

Cap1 = (T12 * T21 - T22 * T11) / (B(2) * T12 - B(1) * T22);
Cap2 = (T12 * T21 - T22 * T11) / (B(1) * T21 - B(2) * T11);

R1 = tau1 / Cap1;
R2 = tau2 / Cap2;

%%
% % alternatively, solve using 
% syms T11t T22t T21t T12t R1t R2t Cap1t Cap2t;
% C1 = C(1);
% C2 = C(2);
% A11 = A(1,1);
% A12 = A(1,2);
% A21 = A(2,1);
% A22 = A(2,2);
% B1 = B(1);
% B2 = B(2);
% eqns = [C1 * T11t + C2 * T21t == 1,...
%         C1 * T12t + C2 * T22t == 1,...
%         B1 == T11t / Cap1t + T12t / Cap2t,...
%         B2 == T21t / Cap1t + T22t / Cap2t,...
%         A11 * T11t + A12 * T21t == -T11t / (R1t * Cap1t),...
%         A11 * T12t + A12 * T22t == -T12t / (R2t * Cap2t),...
%         A21 * T11t + A22 * T21t == -T21t / (R1t * Cap1t),...
%         A21 * T12t + A22 * T22t == -T22t / (R2t * Cap2t)];
% St = solve(eqns);
% 
% T = [round(St.T11t(1),5) round(St.T12t(1),5);...
%     round(St.T21t(1),5) round(St.T22t(1),5)];
% R1 = round(St.R1t(1),5);
% R2 = round(St.R2t(1),5);
% Cap1 = round(St.Cap1t(1),5);
% Cap2 = round(St.Cap2t(1),5);
% R0 = D;

end