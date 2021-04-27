# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 12:42:10 2021

@author: nevetts
"""


import numpy as np
import matplotlib.pyplot as plt
from magpylib.source.magnet import Box,Cylinder
from magpylib import Collection, displaySystem
from magpylib import source
from scipy.optimize import leastsq
import MakeCoils
import matplotlib

def equations(z,a1,b1,c1): 
    # z[0] = thetai,z[1]=xi,a1=thetai-1 , b1=xi-1, c1=L/N

    return(np.cos(a1)-np.cos(z[0])-c1/b1,z[1]-b1-c1/np.cos(a1))

def check(d1,d2,t1,t2,c1):
    if np.abs(np.cos(2*t2)-c1/d2)<1e-8 and np.abs(t2-t1-d2/2-d1/2)<1e-8:
        return True
    else:
        return False
    
def fitFUNC(p,x):
    y_ = p[0]*x+p[1]
    return y_

def timeresiduals(p,x,y):
    """
    returns error between some data and the prediction from the above model
    """
    ymodel=fitFUNC(p,x)
    err=y-ymodel
    return err

matplotlib.rcParams.update({'font.size': 12})

Nt=8

rcoil=133.35/2

print('min Quad gap is ',rcoil/Nt/2)

Hcoil=226*1.5

dt=0.0
delta=(dt)/rcoil/2

Dmin=rcoil/Nt/2#3.3393
Rsample=40
Nphi=50

rot=np.pi/4
c,Coils,Backs = MakeCoils.getQuadCol(rcoil,Hcoil,Dmin,Nt,Nphi,rot)
vertices=c.sources[0].vertices

xs = np.linspace(-10,10,20)
POS = np.array([(x,0,0) for x in xs])
Bs = c.getB(POS)    #<--VECTORIZED

G=np.max(Bs[:,1])# G/A/cm
guess=[-G,0]
plsq = leastsq(timeresiduals, guess, args=(xs,Bs[:,1]),full_output=1) # python finds the No and Rate which minimize the error^2
G=np.abs(plsq[0][0])

vs=np.array(vertices)
pointlength=np.sum((vs[1:]-vs[:-1])**2,axis=1)**0.5
wirelength=np.sum(pointlength)

wirewidth=Dmin-0.1
rho=1.6e-8 # Ohm mm
t=0.035e-3 # mm
nsqu=wirelength/wirewidth
Resistance = rho/t*nsqu
print('Quad R=',Resistance,' Ohms')

Gtarget=0.01
I=Gtarget/G
print('Gradient (G/cm)',G*100)
print('I (A) = ',I)
print('CW Power (W) = ',I**2*Resistance)



fig = plt.figure(figsize=(9,5))
#for i in range(0,npts):
#### Fix bottom!!!
plt.subplot(221)
for i in range(0,len(Coils)):
    MakeCoils.MakeQuadDXF(Coils[i],rcoil,'QC'+str(i))
    plt.plot(rcoil*Coils[i][:,1],Coils[i][:,0],'b-')
    plt.xlim(-20,480)
    plt.ylabel('x (mm)')
    plt.gca().set_title('$G_{zy}$')
    
plt.subplot(224)
for i in range(0,len(Backs)):
    MakeCoils.MakeQuadDXF(Backs[i],rcoil,'QB'+str(i))
    plt.plot(rcoil*Backs[i][:,1],Backs[i][:,0],'r--')
    plt.gca().set_title('excess connections')
    
    
plt.subplot(223)
for i in range(0,len(Coils)):
    Coils[i][:,1]+=rot
    MakeCoils.MakeQuadDXF(Coils[i],rcoil,'Q2C'+str(i))
    plt.plot(rcoil*(Coils[i][:,1]),Coils[i][:,0],'g-')
    plt.xlim(-20,480)
    plt.gca().set_title('$G_{zz}$')
    plt.xlabel('circumference (mm)')
    plt.ylabel('x (mm)')
    
plt.subplot(224)
for i in range(0,len(Backs)):
    Backs[i][:,1]+=rot
    MakeCoils.MakeQuadDXF(Backs[i],rcoil,'Q2B'+str(i))
    plt.plot(rcoil*(Backs[i][:,1]),Backs[i][:,0],'r--')


Nt=20
nc=4000
Nphi=8

wirewidth=3.0
gap=0.1
noff=-0.5
Hmax=Hcoil/2
print('theta=0 spacing, wirewidth Zgrad ',Hmax/Nt,1/Nt,wirewidth)

c,fc,Coils,Backs = MakeCoils.getZGradCol(rcoil,Hmax,Nt,wirewidth,gap,nc)
vertices=c.sources[0].vertices

xs = np.linspace(-10,10,20)
POS = np.array([(0,0,x) for x in xs])
Bs = c.getB(POS)    #<--VECTORIZED

G=np.max(Bs[:,1])# G/A/cm
guess=[-G,0]
plsq = leastsq(timeresiduals, guess, args=(xs,Bs[:,1]),full_output=1) # python finds the No and Rate which minimize the error^2
G=np.abs(plsq[0][0])

vs=np.array(vertices)
pointlength=np.sum((vs[1:]-vs[:-1])**2,axis=1)**0.5
wirelength=np.sum(pointlength)

rho=1.6e-8 # Ohm mm
t=0.035e-3 # mm
nsqu=wirelength/wirewidth
Resistance = rho/t*nsqu
print('Z R=',Resistance,' Ohms')

Gtarget=0.01
I=Gtarget/G
print('Gradient (G/cm)',G*100)
print('I (A) = ',I)
print('CW Power (W) = ',I**2*Resistance)

plt.subplot(222)
for i in range(0,len(Coils)):
    Coils[i][:,1]+=np.pi/2
    MakeCoils.MakeZCoilDXF(Coils[i],rcoil,'ZC'+str(i))
    plt.plot(rcoil*(Coils[i][:,1]),Coils[i][:,0],'k-')
    plt.gca().set_title('$G_{zx}$')
    plt.xlim(-20,480)

plt.subplot(224)
for i in range(0,len(Backs)):
    Backs[i][:,1]+=np.pi/2
    MakeCoils.MakeZCoilDXF(Backs[i],rcoil,'ZB'+str(i))
    plt.plot(rcoil*(Backs[i][:,1]),Backs[i][:,0],'r--')
    plt.xlim(-20,480)



#plt.ylabel('x (mm)')
plt.xlabel('circumference (mm)')
plt.tight_layout()
plt.show()
