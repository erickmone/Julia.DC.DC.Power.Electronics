using ControlSystems
using Plots
using LaTeXStrings

L  = 2.7648e-3;
C  = 1.666667e-6;
R  = 144;
E  = 48;
D  = 0.6;
Dp = 1-D;
X2 = E/Dp;
X1 = E/(R*((Dp^2)));

# Average state-space model, small signal method

mA  = [0 -Dp/L; Dp/C -1/(C*R)];
mB  = [X2/L; -X1/C];
mC1 = [1 0];
mC2 = [0 1];
mD  = 0;

P0= ss(mA,mB,mC2,mD);
P1 = ss(mA,mB,mC1,mD);

# Computation of Transfer Functions

G0 = tf(P0);
G1 = tf(P1);

f1=bodeplot(LTISystem[G0], label=L"G_0");
f2=bodeplot(LTISystem[G1], label=L"G_1");
f3=bodeplot(LTISystem[G0], label=L"G_{vd}",c=:blue);
f4=bodeplot(LTISystem[G1], label=L"G_{id}",c=:orange);

ωLC = Dp/(sqrt(L*C));
Q = R*(sqrt(C/L));
ωz = Q/ωLC;

