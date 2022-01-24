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

tau1
tau2

T11 = 1 / C(1) - (C(2) * A(1,1) * tau1 + C(2)) / ((C(1) * C(2) * A(1,1) - C(1)^2 * A(1,2)) * tau1 + C(1) * C(2));
T12 = 1 / C(1) - (C(2) * A(1,1) * tau2 + C(2)) / ((C(1) * C(2) * A(1,1) - C(1)^2 * A(1,2)) * tau2 + C(1) * C(2));
T21 = (A(1,1) * tau1 + 1) / ((C(2) * A(1,1) - C(1) * A(1,2)) * tau1 + C(2));
T22 = (A(1,1) * tau2 + 1) / ((C(2) * A(1,1) - C(1) * A(1,2)) * tau2 + C(2));
T = [T11 T12; T21 T22];

Cap1 = (T12 * T21 - T22 * T11) / (B(2) * T12 - B(1) * T22);
Cap2 = (T12 * T21 - T22 * T11) / (B(1) * T21 - B(2) * T11);

R1 = tau1 / Cap1;
R2 = tau2 / Cap2;
end