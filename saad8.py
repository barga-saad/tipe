from numpy import *
from matplotlib.pyplot import *

vdc=400
p=1
tp=20
eps=8.85*10**(-12)
cd=10,1
dp=2*10**(-6)
'''
dg=10 #le parametre qu'on veut visualiser son influence sur l'efficacité
'''
a=1.25*10**(-3)
b=150*10**(-3)
l=10
er=10
ld=6.61*10**(-8)*(tp/293)*(101300/p)
Cm=1+2.54*(ld/dp)+(ld/dp)*exp(-0.55*dp/ld)
u=2.4*10**(-5)

def champ_ectrostatique(r=1.25*10**(-4)):
    E=vdc/(r*log(b/a))
    return E


def charge_particule(t):
    qps=pi*eps*(dp**2)*E*((3*er)/(er+2))
    T=(4*eps)/cd
    qp=qps*(t/t+T)
    return qp

def vitesse_de_migration(dg):
    E=champ_ectrostatique(r=1.25*10**(-4))
    we=dp*Cm*(E**2)*((eps*er)/u*(e+2))
    return we


def efficacite(dg):
    E=champ_ectrostatique(r=1.25*10**(-4))
    we=vitesse_de_migration(dg)
    s=2*pi*b*l
    n=1-exp(-we*dp*(s/dg))
    return n
figure(1)
X=linspace(0,60)
Y=efficacite(X)
plot(X,Y,'b*-')
xlabel('dg')
ylabel("l'efficacité")
title("l'efficacité en fonction de debit de gaz")
legend()
grid()
show()