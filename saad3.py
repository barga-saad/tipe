from numpy import *
from matplotlib.pyplot import *

#initialisation de données numerique
vdc=4000
p=10**5
tp=150
eps=8.85*10**(-12)
cd=10,1
dp=2*10**(-6)
dg=10
a=1.25*10**(-3)
'''
b=150*10**(-3) #le parametre qu'on veut visualiser son influence sur l'efficacité
'''
l=10
er=10
ld=6.61*10**(-8)*(tp/293)*(101300/p)
Cm=1+2.54*(ld/dp)+(ld/dp)*exp(-0.55*dp/ld)
u=2.4*10**(-5)

def champ_ectrostatique(r,b):
    E=vdc/(r*log(b/a))
    return E


def charge_particule(t):
    qps=pi*eps*(dp**2)*E*((3*er)/(er+2))
    T=(4*eps)/cd
    qp=qps*(t/t+T)
    return qp

def vitesse_de_migration(E):
    we=dp*Cm*(E**2)*((eps*er)/u*(e+2))
    return we


def efficacite(b,dp,l):
    E=champ_ectrostatique(b,0.003)
    we=vitesse_de_migration(E)
    s=2*pi*b*l
    n=1-exp(-we*dp*(s/dg))
    return n
    

figure(3)
x=linspace(0.5,2)
y1=efficacite(0.12,x,50)
y2=efficacite(1.2,x,50)
y3=efficacite(12,x,50)
plot(x,y1,'*-',label='B=0.12m')
plot(x,y2,'+-',label='B=1.2m')
plot(x,y3,'--',label='BL=12m')
ylabel('efficacité')
xlabel('diametre de la particule')
title("l'efficacité en fonction de dp")
legend()
grid()
show()