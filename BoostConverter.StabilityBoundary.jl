using ControlSystems
using Plots
using Polynomials 

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

Pid = Gid;
doplot = true
form = :parallel
kp, ki, f1 = stabregionPID(Pid,range(0,stop=20e3,length=9000); doplot, form); 
