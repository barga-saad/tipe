from numpy import *
from matplotlib.pyplot import *

#initialisation de données numerique
vdc=400
p=10**5
tp=150
eps=8.85*10**(-12)
cd=10,1
dp=2*10**(-6)
dg=10
a=1.25*10**(-3)
b=150*10**(-3)
l=10
er=10
ld=6.61*10**(-8)*(tp/293)*(101300/p)
Cm=1+2.54*(ld/dp)+(ld/dp)*exp(-0.55*dp/ld)
u=2.4*10**(-5)

def champ_ectrostatique(r):
    E=vdc/(r*log(b/a))
    return E


def charge_particule(t):
    qps=pi*eps*(dp**2)*E*((3*er)/(er+2))
    T=(4*eps)/cd
    qp=qps*(t/t+T)
    return qp

def vitesse_de_migration(r):
    E=champ_ectrostatique(r)
    we=dp*Cm*(E**2)*((eps*er)/u*(e+2))
    return we


def efficacite(r):
    E=champ_ectrostatique(r)
    we=vitesse_de_migration(r)
    s=2*pi*b*l
    n=1-exp(-we*dp*(s/dg))
    return n
    
X=linspace(0.001,0.12)
Y=champ_ectrostatique(X)
plot(X,Y,'b')
xlabel=('champ electrique')
ylabel=('r')
title('le champ electrique en fonction de r')
legend()
grid()
show()