function [Re0, Re1, Re2, Re3, Ca1, Ca2, Ca3, T] = extract_3rc_params(A,B,C,D)
% return the resistor and capacitor values in a 3RC equivalent circuit
% battery model. Also return the transformation matrix T that converts
% between the obtained state space and the desired 3RC state space

% check that all matrices are of the correct size
assert(all(size(A)==[3 3]))
assert(all(size(B)==[3 1]))
assert(all(size(C)==[1 3]))
assert(all(size(D)==[1 1]))

% perform parameter extraction using MATLAB's equation solver
syms T11 T12 T13 T21 T22 T23 T31 T32 T33 R1 R2 R3 Cap1 Cap2 Cap3;
C1 = C(1);
C2 = C(2);
C3 = C(3);
A11 = A(1,1);
A12 = A(1,2);
A13 = A(1,3);
A21 = A(2,1);
A22 = A(2,2);
A23 = A(2,3);
A31 = A(3,1);
A32 = A(3,2);
A33 = A(3,3);
B1 = B(1);
B2 = B(2);
B3 = B(3);
eqns = [C1 * T11 + C2 * T21 + C3 * T31 == 1,...
        C1 * T12 + C2 * T22 + C3 * T32 == 1,...
        C1 * T13 + C2 * T23 + C3 * T33 == 1,...
        B1 == T11 / Cap1 + T12 / Cap2 + T13 / Cap3,...
        B2 == T21 / Cap1 + T22 / Cap2 + T23 / Cap3,...
        B3 == T31 / Cap1 + T32 / Cap2 + T33 / Cap3,...
        A11 * T11 + A12 * T21 + A13 * T31 == -T11 / (R1 * Cap1),...
        A21 * T11 + A22 * T21 + A23 * T31 == -T21 / (R1 * Cap1),...
        A31 * T11 + A32 * T21 + A33 * T31 == -T31 / (R1 * Cap1),...
        A11 * T12 + A12 * T22 + A13 * T32 == -T12 / (R2 * Cap2),...
        A21 * T12 + A22 * T22 + A23 * T32 == -T22 / (R2 * Cap2),...
        A31 * T12 + A32 * T22 + A33 * T32 == -T32 / (R2 * Cap2),...
        A11 * T13 + A12 * T23 + A13 * T33 == -T13 / (R3 * Cap3),...
        A21 * T13 + A22 * T23 + A23 * T33 == -T23 / (R3 * Cap3),...
        A31 * T13 + A32 * T23 + A33 * T33 == -T33 / (R3 * Cap3)];
St = solve(eqns);

T = [round(St.T11(1),5) round(St.T12(1),5) round(St.T13(1),5);...
    round(St.T21(1),5) round(St.T22(1),5) round(St.T23(1),5);...
    round(St.T31(1),5) round(St.T32(1),5) round(St.T33(1),5)];

Re0 = D;
Re1 = round(St.R1(1),5);
Re2 = round(St.R2(1),5);
Re3 = round(St.R3(1),5);

Ca1 = round(St.Cap1(1),5);
Ca2 = round(St.Cap2(1),5);
Ca3 = round(St.Cap3(1),5);


end