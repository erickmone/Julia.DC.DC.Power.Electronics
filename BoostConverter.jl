using ControlSystems
using Plots

L  = 2.7648e-3;
C  = 1.666667e-6;
R  = 144;
E  = 48;
D  = 0.6;
Dp = 1-D;
X1 = E/Dp;
X2 = E/(R*((Dp^2)));

# Average state-space model, small signal method

mA  = [0 -Dp/L; Dp/C -1/(C*R)];
mB  = [X2/L; -X1/C];
mC1 = [1 0];
mC2 = [0 1];
mD  = 0;

P0= ss(mA,mB,mC1,mD);
P1 = ss(mA,mB,mC2,mD);

# Computation of Transfer Functions

G0 = tf(P0);
G1 = tf(P1);

bodeplot(LTISystem[G0])
