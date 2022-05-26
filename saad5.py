from numpy import *
from matplotlib.pyplot import *

vdc=400*10**4
p=1
'''
tp=150    #le parametre qu'on veut visualiser son influence sur l'efficacité
'''
eps=8.85*10**(-12)
cd=10,1
dp=2*10**(-6)
dg=10
a=1.25*10**(-3)
b=150*10**(-3)
l=10
er=10
'''
ld=6.61*10**(-8)*(tp/293)*(101300/p)     #VARIABLE
Cm=1+2.54*(ld/dp)+(ld/dp)*exp(-0.55*dp/ld) #VARIABLE
'''
u=2.4*10**(-5)

def champ_ectrostatique(r):
    E=vdc/(r*log(b/a))
    return E


def charge_particule(t):
    qps=pi*eps*(dp**2)*E*((3*er)/(er+2))
    T=(4*eps)/cd
    qp=qps*(t/t+T)
    return qp

def vitesse_de_migration(tp):
    E=400*10**4
    ld=6.61*10**(-8)*(tp/293)*(101300/p)
    Cm=1+2.54*(ld/dp)+(ld/dp)*exp(-0.55*dp/ld)
    we=dp*Cm*(E**2)*((eps*er)/u*(e+2))
    return we


def efficacite(tp):
    E=400*10**4
    we=vitesse_de_migration(tp)
    s=2*pi*b*l
    n=1-exp(-we*dp*(s/dg))
    return n
figure(1)
X=linspace(20,50)
Y=efficacite(X)
plot(X,Y,'b')
xlabel('temperature')
ylabel("l'efficacité")
title("l'efficacité en fonction de la temperature")
legend()
grid()
show()