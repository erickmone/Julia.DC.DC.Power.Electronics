using ControlSystems
using Plots
using LaTeXStrings
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

n1 = numI[1];
n0 = numI[2];
d2 = den[1];
d1 = den[2];
d0 = den[3];

PnumI = Polynomial([n0,n1], :s);
ZnumI = roots(PnumI);

Pden = Polynomial([d0,d1,d2], :s);
Rden = roots(Pden);

rlocusplot(Gvd,xlims=(-4e4,4e4),ylims=(-2.5e4,2.5e4),)
rlocusplot!(Gid,xlims=(-4e4,4e4),ylims=(-2.5e4,2.5e4),)

σ1 = 1000;
Σnum1 = [n1, n0-(n1*σ1)];
Σden1 = [d2, d1-(2*d2*σ1), (σ1^2)-(d1*σ1)+d0];
Gσ1 = tf(Σnum1,Σden1);

σ2 = 2000;
Σnum2 = [n1, n0-(n1*σ2)];
Σden2 = [d2, d1-(2*d2*σ2), (σ2^2)-(d1*σ2)+d0];
Gσ2 = tf(Σnum2,Σden2);

σ3 = 3000;
Σnum3 = [n1, n0-(n1*σ3)];
Σden3 = [d2, d1-(2*d2*σ3), (σ3^2)-(d1*σ3)+d0];
Gσ3 = tf(Σnum3,Σden3);

σ4 = 4000;
Σnum4 = [n1, n0-(n1*σ4)];
Σden4 = [d2, d1-(2*d2*σ4), (σ4^2)-(d1*σ4)+d0];
Gσ4 = tf(Σnum4,Σden4);

Pid = Gid;
Pσ1 = Gσ1;
Pσ2 = Gσ2;
Pσ3 = Gσ3;
Pσ4 = Gσ4;


doplot = true
form = :parallel

kp, ki, f1 = stabregionPID(Pid,range(0,stop=20e3,length=9000); doplot, form); 
kp1, ki1, f2 = stabregionPID(Pσ1,range(0,stop=20e3,length=9000); doplot, form); 
kp2, ki2, f3 = stabregionPID(Pσ2,range(0,stop=20e3,length=9000); doplot, form); 
kp3, ki3, f4 = stabregionPID(Pσ3,range(0,stop=20e3,length=9000); doplot, form); 
kp4, ki4, f5 = stabregionPID(Pσ4,range(0,stop=20e3,length=9000); doplot, form); 

plot(kp,ki)
plot!(kp1,ki1)
plot!(kp2,ki2)
plot!(kp3,ki3)
plot!(kp4,ki4,xlims=(-0.1, 0.1), ylims=(0,1500))


# https://goropikari.github.io/PlotsGallery.jl/src/parametric2d.html