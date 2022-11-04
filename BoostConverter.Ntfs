using ControlSystems
using Plots
using LaTeXStrings

L  = 2.7648e-3;
C  = 1.666667e-6;
R  = 144;
E  = 48;
D  = 0.6;
Dp = 1-D;

V = E/Dp; # Equilibrium of output voltage

# Construction of the Transfer Functions

Gd1 = (2*V)/((Dp^2)*R);
Gd0 = V/Dp;

ωy = 2/(R*C);
ωz = -((Dp^2)*R)/L;
ωo = Dp/(sqrt(L*C));
Q  = Dp*R*sqrt(C/L);

numI = [Gd1*((ωo^2)/ωy), Gd1*(ωo^2)];
numV = [Gd0*((ωo^2)/ωz), Gd0*(ωo^2)];
den = [1, ωo/Q, ωo^2];

Gid = tf(numI,den); # TF inductor current-to-duty-cycle
Gvd = tf(numV,den); # TF output voltage-to-duty-cycle


n1 = numI[1];
n0 = numI[2];
d2 = den[1];
d1 = den[2];
d0 = den[3];

rlocusplot(Gvd,xlims=(-4e4,4e4),ylims=(-2.5e4,2.5e4),)
rlocusplot!(Gid,xlims=(-4e4,4e4),ylims=(-2.5e4,2.5e4),)
