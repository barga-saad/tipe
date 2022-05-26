from numpy import *
from matplotlib.pyplot import *

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

def vitesse_de_migration(E):
    we=dp*Cm*(E**2)*((eps*er)/u*(e+2))
    return we


def efficacite(E,dp,l):
    we=vitesse_de_migration(E)
    s=2*pi*b*l
    n=1-exp(-we*dp*(s/dg))
    return n
    
    
figure(1)
X=linspace(0.001,0.12)
Y=champ_ectrostatique(X)
plot(X,Y,'b')
ylabel('champ electrique')
xlabel('r')
title('le champ electrique en fonction de r')
legend()
grid()
show()

figure(2)
x=linspace(0.5,2)
y1=efficacite(40000,x,10)
y2=efficacite(40000,x,20)
y3=efficacite(40000,x,30)
y4=efficacite(40000,x,50)
plot(x,y1,label='L=10m')
plot(x,y2,'*-',label='L=20m')
plot(x,y3,'+-',label='L=30m')
plot(x,y4,'--',label='L=50m')
ylabel('efficacité')
xlabel('diametre de la particule')
title("l'efficacité en fonction de dp")
legend()
grid()
show()
