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
import ezdxf
from ezdxf import units
from ezdxf.addons import r12writer


def getIwire(rc,Hm,x_,t_,ww,wg,turn,tlast=np.pi/2):
    
    Dm=ww+wg
    if turn ==0:
        ds=2*rc*(tlast-t_[:,turn])
        I2=np.max(np.where(ds>Dm))
    else:
        ds=rc*(tlast-t_[:,turn])
        I2=np.max(np.where(ds>Dm))
    
    I1 = np.max(np.where(x_[:,turn]<Hm-turn*Dm))
    
    If=np.min([I1,I2])
    
    return If

def getZGradCol(rcoil,Hmax,Nt,wirewidth,gap,nc=1500,Nphi=8,noff=-0.5):
    
    flatcoords=[]
    backcoords=[]
    Q1=[]
    Q2=[]
    Q3=[]
    Q4=[]
    
    B1=[]
    B2=[]
    B3=[]
    B4=[]

    Dmin=wirewidth+gap
    Hcoil=Hmax-(Nt+noff+1)*Dmin
    
    x0=Hcoil/Nt
    theta0=0
    X=np.zeros((nc,Nt))
    T=np.zeros((nc,Nt))
    Xw=np.zeros(nc)
    
    for j in range(0,Nt):
        X[0,j]=(j+1)*x0+noff*x0
        
        
        for i in range(1,nc):
          
            X[i,j]=X[i-1,j]+Hcoil/nc/np.cos(T[i-1,j])
            T[i,j]=np.arccos(np.cos(T[i-1,j])-Hcoil/nc/X[i,j])
            
            
    vertices=[]
    Tedge=np.pi/2
    xedge=Hmax
    
    ds=2*rcoil*(Tedge-T[:,0])
    ITmaxnext=getIwire(rcoil,Hmax,X,T,wirewidth,gap,0)
    
    
    ### Q1 
    for j in range(0,Nt):

        dsi=rcoil*(Tedge-T[:,j])
        dx=Hmax-(j-1)*Dmin-X[ITmaxnext,j]    
        
        ITmax=ITmaxnext
        #print('azimuthal,axial sep ',j,dsi[ITmax],dx)
        Tedge=T[ITmax,j]
        xedge=X[ITmax,j]
        
        xi=rcoil*np.cos(T[ITmax,j])
        yi=rcoil*np.sin(T[ITmax,j])
        zi=Hmax-j*Dmin
        vertices.append([xi,yi,zi])
        flatcoords.append([zi,T[ITmax,j]])
        Q1.append([zi,T[ITmax,j]])
        
        
        for i in range(ITmax-1,0,-1):
            xi=rcoil*np.cos(T[i,j])
            yi=rcoil*np.sin(T[i,j])
            zi=X[i,j]
        
            vertices.append([xi,yi,zi])
            flatcoords.append([zi,T[i,j]])
            Q1.append([zi,T[i,j]])
            
        for i in range(0,ITmax-1):
            xi=rcoil*np.cos(-1*T[i,j])
            yi=rcoil*np.sin(-1*T[i,j])
            zi=X[i,j]

            vertices.append([xi,yi,zi])
            flatcoords.append([zi,-T[i,j]])
            Q1.append([zi,-T[i,j]])
        
        xi=rcoil*np.cos(-T[ITmax,j])
        yi=rcoil*np.sin(-T[ITmax,j])
        zi=Hmax-j*Dmin
        vertices.append([xi,yi,zi])
        flatcoords.append([zi,-T[ITmax,j]])
        Q1.append([zi,-T[ITmax,j]])
        
        if j!=Nt-1:
            
            jnext=j+1
            ITmaxnext=getIwire(rcoil,Hmax,X,T,wirewidth,gap,jnext,Tedge)
            phi=np.linspace(-T[ITmax,j],T[ITmaxnext,jnext],Nphi)
            # around circumference
            for k in range(0,Nphi):
                vertices.append([rcoil*np.cos(phi[k]),rcoil*np.sin(phi[k]),zi])
                flatcoords.append([zi,phi[k]])
                Q1.append([zi,phi[k]])
        else:
            jnext=j
            zi=Hmax
        

    ###### Q1 - Q2
    Istart = getIwire(rcoil,Hmax,X,T,wirewidth,gap,0)#np.max(np.where(X[:,0]<Hmax)) 
    phi=np.linspace(-T[ITmax,Nt-1],np.pi-T[Istart,0],Nphi)
    
    for k in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[k]),rcoil*np.sin(phi[k]),zi])
    
    flatcoords.append([zi,phi[0]])
    flatcoords.append([zi,phi[-1]])
    
    B1.append([X[ITmax,Nt-1],-T[ITmax,Nt-1]])
    B1.append([Hmax,phi[0]])
    B1.append([Hmax,phi[-1]])
    
    ##Q2
    ITmaxnext=getIwire(rcoil,Hmax,X,T,wirewidth,gap,0)
    for j in range(0,Nt):
        
        ITmax=ITmaxnext
        Tedge=T[ITmax,j]
        xedge=X[ITmax,j]
        
        xi=rcoil*np.cos(np.pi-T[ITmax,j])
        yi=rcoil*np.sin(np.pi-T[ITmax,j])
        zi=(Hmax-j*Dmin)
        vertices.append([xi,yi,zi])
        flatcoords.append([zi,np.pi-T[ITmax,j]])
        Q2.append([zi,np.pi-T[ITmax,j]])
        
        
        for i in range(ITmax-1,0,-1):
            xi=rcoil*np.cos(np.pi-T[i,j])
            yi=rcoil*np.sin(np.pi-T[i,j])
            zi=X[i,j]
                
            vertices.append([xi,yi,zi])
            flatcoords.append([zi,np.pi-T[i,j]])
            Q2.append([zi,np.pi-T[i,j]])
            
        for i in range(0,ITmax-1):
            xi=rcoil*np.cos(np.pi+T[i,j])
            yi=rcoil*np.sin(np.pi+T[i,j])
            zi=X[i,j]
        
        
            vertices.append([xi,yi,zi])
            flatcoords.append([zi,np.pi+T[i,j]])
            Q2.append([zi,np.pi+T[i,j]])
            
        xi=rcoil*np.cos(np.pi+T[ITmax,j])
        yi=rcoil*np.sin(np.pi+T[ITmax,j])
        zi=(Hmax-j*Dmin)
        vertices.append([xi,yi,zi])
        flatcoords.append([zi,np.pi+T[ITmax,j]])
        Q2.append([zi,np.pi+T[ITmax,j]])
        
        if j!=Nt-1:
            
            jnext=j+1
            ITmaxnext=getIwire(rcoil,Hmax,X,T,wirewidth,gap,jnext,Tedge)
            phi=np.linspace(np.pi+T[ITmax,j],np.pi-T[ITmaxnext,jnext],Nphi)
        
            # around circumference
            for k in range(0,Nphi):
                vertices.append([rcoil*np.cos(phi[k]),rcoil*np.sin(phi[k]),zi])
                flatcoords.append([zi,phi[k]])
                Q2.append([zi,phi[k]])
        else:
            jnext=j
            zi=Hmax
    
    ##Q2 - Q3
    backoff=2*Dmin
    backphi=backoff/rcoil
    phi=np.linspace(np.pi+T[ITmax,Nt-1],np.pi-T[Istart,0],Nphi)
    
    for k in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[k]),rcoil*np.sin(phi[k]),zi])
            
    flatcoords.append([zi,phi[0]])
    flatcoords.append([zi,phi[-1]])
    
    B2.append([X[ITmax,Nt-1],phi[0]])# end of Q2
    B2.append([Hmax,phi[0]])# back to top of coil
    B2.append([Hmax,phi[-1]+backphi]) # back to coil almost center
    B2.append([Hmax-backoff,phi[-1]+backphi]) # down a bit
    B2.append([Hmax-backoff,phi[-1]]) # coil center
    B2.append([-Hmax+backoff/2,phi[-1]]) # almost Q3
    
    ## Q3
    ITmaxnext=getIwire(rcoil,Hmax,X,T,wirewidth,gap,0)
    for j in range(0,Nt):        
        ITmax=ITmaxnext
        Tedge=T[ITmax,j]
        
        xi=rcoil*np.cos(np.pi-T[ITmax,j])
        yi=rcoil*np.sin(np.pi-T[ITmax,j])
        zi=-1*(Hmax-j*Dmin)
        vertices.append([xi,yi,zi])
        flatcoords.append([zi,np.pi-T[ITmax,j]])
        Q3.append([zi,np.pi-T[ITmax,j]])
        
        for i in range(ITmax-1,0,-1):
            xi=rcoil*np.cos(np.pi-T[i,j])
            yi=rcoil*np.sin(np.pi-T[i,j])
            zi=-X[i,j]
        
            vertices.append([xi,yi,zi])
            flatcoords.append([zi,np.pi-T[i,j]])
            Q3.append([zi,np.pi-T[i,j]])
            
        for i in range(0,ITmax-1):
            xi=rcoil*np.cos(np.pi+T[i,j])
            yi=rcoil*np.sin(np.pi+T[i,j])
            zi=-X[i,j]
        
            vertices.append([xi,yi,zi])
            flatcoords.append([zi,np.pi+T[i,j]])
            Q3.append([zi,np.pi+T[i,j]])
            
        xi=rcoil*np.cos(np.pi+T[ITmax,j])
        yi=rcoil*np.sin(np.pi+T[ITmax,j])
        zi=-1*(Hmax-j*Dmin)
        vertices.append([xi,yi,zi])
        flatcoords.append([zi,np.pi+T[ITmax,j]])
        Q3.append([zi,np.pi+T[ITmax,j]])
        
        if j!=Nt-1:
            jnext=j+1
            
            ITmaxnext=getIwire(rcoil,Hmax,X,T,wirewidth,gap,jnext,Tedge)
            phi=np.linspace(np.pi+T[ITmax,j],np.pi-T[ITmaxnext,jnext],Nphi)
        
            # around circumference
            for k in range(0,Nphi):
                vertices.append([rcoil*np.cos(phi[k]),rcoil*np.sin(phi[k]),zi])
                flatcoords.append([zi,phi[k]])
                Q3.append([zi,phi[k]])
        else:
            jnext=j
            zi=-Hmax
    
    
    ##Q3-Q4
    phi=np.linspace(np.pi+T[ITmax,Nt-1],T[Istart,0],Nphi)
    for k in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[k]),rcoil*np.sin(phi[k]),zi])
            
    flatcoords.append([zi,phi[0]])
    flatcoords.append([zi,phi[-1]])
    
    B3.append([-X[ITmax,Nt-1],phi[0]])# end of Q3
    B3.append([-Hmax,phi[0]])# back to bottm of coil
    B3.append([-Hmax,phi[-1]]) # back Q4
    
    
    ##Q4
    ITmaxnext=getIwire(rcoil,Hmax,X,T,wirewidth,gap,0)
    for j in range(0,Nt):

        ITmax=ITmaxnext
        Tedge=T[ITmax,j]
        xedge=X[ITmax,j]
        
        xi=rcoil*np.cos(T[ITmax,j])
        yi=rcoil*np.sin(T[ITmax,j])
        zi=-1*(Hmax-j*Dmin)
        vertices.append([xi,yi,zi])
        flatcoords.append([zi,T[ITmax,j]])
        Q4.append([zi,T[ITmax,j]])
        
        
        for i in range(ITmax-1,0,-1):
            xi=rcoil*np.cos(T[i,j])
            yi=rcoil*np.sin(T[i,j])
            zi=-X[i,j]
        
            vertices.append([xi,yi,zi])
            flatcoords.append([zi,T[i,j]])
            Q4.append([zi,T[i,j]])
            
        for i in range(0,ITmax-1):
            xi=rcoil*np.cos(-1*T[i,j])
            yi=rcoil*np.sin(-1*T[i,j])
            zi=-X[i,j]
        
            vertices.append([xi,yi,zi])
            flatcoords.append([zi,-T[i,j]])
            Q4.append([zi,-T[i,j]])
            
        xi=rcoil*np.cos(-T[ITmax,j])
        yi=rcoil*np.sin(-T[ITmax,j])
        zi=-1*(Hmax-j*Dmin)
        vertices.append([xi,yi,zi])
        flatcoords.append([zi,-T[ITmax,j]])
        Q4.append([zi,-T[ITmax,j]])
        
        if j!=Nt-1:
            jnext=j+1
            
            ITmaxnext=getIwire(rcoil,Hmax,X,T,wirewidth,gap,jnext,Tedge)
            phi=np.linspace(-T[ITmax,j],T[ITmaxnext,jnext],Nphi)
            # around circumference
            for k in range(0,Nphi):
                vertices.append([rcoil*np.cos(phi[k]),rcoil*np.sin(phi[k]),zi])
                flatcoords.append([zi,phi[k]])
                Q4.append([zi,phi[k]])
        else:
            jnext=j
            zi=-Hmax
    
    ##### Q4 - Q1
    
    phi=np.linspace(-T[ITmax,Nt-1],T[Istart,0],Nphi)
    for k in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[k]),rcoil*np.sin(phi[k]),zi])
            
    flatcoords.append([zi,phi[0]])
    flatcoords.append([zi,phi[-1]])
            
    vertices.append([rcoil*np.cos(T[Istart,0]),rcoil*np.sin(T[Istart,0]),Hmax])
    flatcoords.append([Hmax,T[Istart,0]])
    
    B4.append([-X[ITmax,Nt-1],phi[0]])# end of Q4
    B4.append([-Hmax,phi[0]])# back to bottom of coil
    B4.append([-Hmax,phi[-1]-backphi/2]) # back to coil almost center
    B4.append([-Hmax+backoff/2,phi[-1]-backphi/2]) # up a bit
    B4.append([-Hmax+backoff/2,phi[-1]]) # over to coil center
    B4.append([Hmax-backoff/2,phi[-1]]) # up the side!
    B4.append([Hmax-backoff/2,phi[-1]+backphi]) # over a bit
    B4.append([Hmax,phi[-1]+backphi]) # and finally out
        
    
    B1=np.array(B1)
    B2=np.array(B2)
    B3=np.array(B3)
    B4=np.array(B4)
    
    Q1=np.array(Q1)
    Q2=np.array(Q2)
    Q3=np.array(Q3)
    Q4=np.array(Q4)
    
    flatcoords=np.array(flatcoords)
    Backside=[B1,B2,B3,B4]
    Coils=[Q1,Q2,Q3,Q4]
    cd = source.current.Line(curr=1,vertices=vertices)
    # create collection
    c = Collection(cd)
    c.rotate(90,(0,0,1))
    
    return [c,flatcoords,Coils,Backside]


def getQuadCol(rcoil,Hcoil,Dmin,Nt,Nphi=50,grotate=0):


    
    c0=2*Nt#/(1-np.sin(2*delta))
    b0=0
    betas=np.zeros(Nt)
    wb=np.zeros(Nt)
    
    
    for i in range(0,Nt):
    
        bp = np.arcsin(1/Nt + np.sin(2*b0)) /2
        betas[i]=bp
        wi=(bp-b0)/2+b0
        wb[i]=wi
        b0=bp
    
    
    brev=np.pi/2-betas
    betas=np.append(betas,brev)
    betas=np.hstack([betas,betas+np.pi/2,betas+np.pi,betas+3*np.pi/2])
    
    wb=wb
    wrev1=np.pi/2-wb
    wrev=wrev1[::-1]
    wb=np.append(wb,wrev)
    #wb=np.hstack([wb,wb+np.pi/2,wb+np.pi,wb+3*np.pi/2])
    
    xwires=rcoil*wb
    pwires=betas*rcoil
    
    vertices=[]
    
    Q1=[]
    Q2=[]
    Q3=[]
    Q4=[]
    Q5=[]
    
    B1=[]
    B2=[]
    B3=[]
    B4=[]
    
    bt=2
    x1=rcoil*np.cos(wb[0]) 
    y1=rcoil*np.sin(wb[0])
    vertices.append([x1,y1,-Hcoil/2])
    Q1.append([-Hcoil/2,wb[0]])
    
    #Q1
    for i in range(0,Nt):
        x1=rcoil*np.cos(wb[i]) 
        y1=rcoil*np.sin(wb[i])
        zb=-Hcoil/2+(i+bt)*Dmin
        zt=Hcoil/2-(i)*Dmin
        
        vertices.append([x1,y1,zb])
        vertices.append([x1,y1,zt]) # up the side!
        
        Q1.append([zb,wb[i]])
        Q1.append([zt,wb[i]]) # up the side!     
        
        phi=np.linspace(wb[i],wrev1[i],Nphi)
        # around circumference
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zt])
            Q1.append([Hcoil/2-i*Dmin,phi[j]])
        
        x2=rcoil*np.cos(wrev1[i])
        y2=rcoil*np.sin(wrev1[i])
        zbnext=-Hcoil/2+(i+bt+1)*Dmin
        
        vertices.append([x2,y2,zbnext]) # down the side!
        Q1.append([zbnext,wrev1[i]]) # down the side!
        
        phi=np.linspace(wrev1[i],wb[i+1],Nphi)
        zbnext=-Hcoil/2+(i+bt+1)*Dmin
        
        
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zbnext])
            Q1.append([zbnext,phi[j]])
    
    ## Q1-Q2
    offset=np.pi/2
    #out1=np.pi/4
    #vertices.append([rcoil*np.cos(out1),rcoil*np.sin(out1),zbnext])
    #Q1.append([zbnext,out1])
    
    zb=-Hcoil/2+(bt)*Dmin
    B1.append([zbnext,phi[j]])
    B1.append([zb,phi[j]])
    phi=np.linspace(phi[j],wrev1[0]+offset,Nphi)
    
    for j in range(0,Nphi):
        vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zb])  
        Q2.append([zb,phi[j]])  
        
    #Q2
    for i in range(0,Nt):
        x1=rcoil*np.cos(wrev1[i]+offset) 
        y1=rcoil*np.sin(wrev1[i]+offset)
        zb=-Hcoil/2+(i+bt)*Dmin
        zt=Hcoil/2-(i)*Dmin
        
        vertices.append([x1,y1,zb])
        vertices.append([x1,y1,zt]) # up the side!
        Q2.append([zb,wrev1[i]+offset])
        Q2.append([zt,wrev1[i]+offset]) # up the side!
        
        phi=np.linspace(wrev1[i]+offset,wb[i]+offset,Nphi)
        # around circumference
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zt])
            Q2.append([zt,phi[j]])
        
        x2=rcoil*np.cos(wb[i]+offset)
        y2=rcoil*np.sin(wb[i]+offset)
        zbnext=-Hcoil/2+(i+bt+1)*Dmin
        
        vertices.append([x2,y2,zbnext]) # down the side!
        Q2.append([zbnext,wb[i]+offset]) # down the side!
        
        phi=np.linspace(wb[i]+offset,wb[::-1][i+1]+offset,Nphi)

        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zbnext])
            Q2.append([zbnext,phi[j]])
    
    
    ## Q2-Q3      
    zb=-Hcoil/2+(bt-1)*Dmin
    
    #vertices.append([rcoil*np.cos(out1+offset),rcoil*np.sin(out1+offset),zbnext])
    #Q2.append([zbnext,out1+offset])
    
    B2.append([zbnext,phi[j]])
    B2.append([zb,phi[j]])
    offset2=np.pi
    phi=np.linspace(phi[j],wb[0]+offset2,Nphi)
    
    for j in range(0,Nphi):
        vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zb])          
        Q3.append([zb,phi[j]])   
        
        
    #Q3    
    for i in range(0,Nt):
        x1=rcoil*np.cos(wb[i]+offset2) 
        y1=rcoil*np.sin(wb[i]+offset2)
        zb=-Hcoil/2+(i+bt)*Dmin
        zt=Hcoil/2-(i)*Dmin
        
        vertices.append([x1,y1,zb])
        vertices.append([x1,y1,zt]) # up the side!
        Q3.append([zb,wb[i]+offset2])
        Q3.append([zt,wb[i]+offset2]) # up the side!   
     
        phi=np.linspace(wb[i]+offset2,wrev1[i]+offset2,Nphi)
        # around circumference
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zt])
            Q3.append([zt,phi[j]])
        
        x2=rcoil*np.cos(wrev1[i]+offset2)
        y2=rcoil*np.sin(wrev1[i]+offset2)
        zbnext=-Hcoil/2+(i+bt+1)*Dmin
        
        vertices.append([x2,y2,zbnext]) # down the side!
        Q3.append([zbnext,wrev1[i]+offset2]) # down the side!
        
        phi=np.linspace(wrev1[i]+offset2,wb[i+1]+offset2,Nphi)
        
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zbnext])
            Q3.append([zbnext,phi[j]])
            
            

    ## Q3-Q4   
    zb=-Hcoil/2+(bt-1)*Dmin
    
    B3.append([zbnext,phi[j]])
    B3.append([zb,phi[j]])
    offset3=3*np.pi/2
    phi=np.linspace(phi[j],wrev1[0]+offset3,Nphi)
    
    for j in range(0,Nphi):
        vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zb])   
        Q4.append([zb,phi[j]])   
    
    #Q4
    for i in range(0,Nt):
        x1=rcoil*np.cos(wrev1[i]+offset3) 
        y1=rcoil*np.sin(wrev1[i]+offset3)
        zb=-Hcoil/2+(i+bt)*Dmin
        zt=Hcoil/2-(i)*Dmin
        
        vertices.append([x1,y1,zb])
        vertices.append([x1,y1,zt]) # up the side!
        Q4.append([zb,wrev1[i]+offset3])
        Q4.append([zt,wrev1[i]+offset3]) # up the side!
        
        phi=np.linspace(wrev1[i]+offset3,wb[i]+offset3,Nphi)
        # around circumference
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zt])
            Q4.append([zt,phi[j]])
        
        x2=rcoil*np.cos(wb[i]+offset3)
        y2=rcoil*np.sin(wb[i]+offset3)
        zbnext=-Hcoil/2+(i+bt+1)*Dmin
        vertices.append([x2,y2,zbnext]) # down the side!
        Q4.append([zbnext,wb[i]+offset3]) # down the side!
        
        phi=np.linspace(wb[i]+offset3,wb[::-1][i+1]+offset3,Nphi)
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zbnext])
            Q4.append([zbnext,phi[j]])
    
    ## Q4-Q1   
    zb=-Hcoil/2
   
    B4.append([zbnext,phi[j]])
    B4.append([zb,phi[j]])
    
    phi=np.linspace(phi[j],wb[1],Nphi)
    for j in range(0,Nphi):
        vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zb])   
        Q5.append([zb,phi[j]])   
    
    spacing=xwires[1:]-xwires[:-1]
    print('argmin ',np.argmin(spacing))
    print('min wire spacing',np.min(np.abs(spacing)))
    print('min wire spacing',2*xwires[0])
    
    B1=np.array(B1)
    B2=np.array(B2)
    B3=np.array(B3)
    B4=np.array(B4)
    
    Q1=np.array(Q1)
    Q2=np.array(Q2)
    Q3=np.array(Q3)
    Q4=np.array(Q4)
    Q5=np.array(Q5)
    

    Backside=[B1,B2,B3,B4]
    Coils=[Q1,Q2,Q3,Q4,Q5]
    
    cd = source.current.Line(curr=1,vertices=vertices)
    
    # create collection
    c = Collection(cd)
    
    return [c,Coils,Backside]


def get2SideQuadCol(rcoil,Hcoil,Dmin,Nt,Nphi=50,grotate=0):


    
    c0=2*Nt#/(1-np.sin(2*delta))
    b0=0
    betas=np.zeros(Nt)
    wb=np.zeros(Nt)
    
    
    for i in range(0,Nt):
    
        bp = np.arcsin(1/Nt + np.sin(2*b0)) /2
        betas[i]=bp
        wi=(bp-b0)/2+b0
        wb[i]=wi
        b0=bp
    
    
    brev=np.pi/2-betas
    betas=np.append(betas,brev)
    betas=np.hstack([betas,betas+np.pi/2,betas+np.pi,betas+3*np.pi/2])
    
    wb=wb
    wrev1=np.pi/2-wb
    wrev=wrev1[::-1]
    wb=np.append(wb,wrev)
    #wb=np.hstack([wb,wb+np.pi/2,wb+np.pi,wb+3*np.pi/2])
    
    xwires=rcoil*wb
    pwires=betas*rcoil
    
    vertices=[]
    
    Q1=[]
    Q2=[]
    Q3=[]
    Q4=[]
    Q5=[]
    
    B1=[]
    B2=[]
    B3=[]
    B4=[]
    
    bt=2
    x1=rcoil*np.cos(wb[0]) 
    y1=rcoil*np.sin(wb[0])
    vertices.append([x1,y1,-Hcoil/2])
    Q1.append([-Hcoil/2,wb[0]])
    
    #Q1
    for i in range(0,Nt):
        x1=rcoil*np.cos(wb[i]) 
        y1=rcoil*np.sin(wb[i])
        zb=-Hcoil/2+(i+bt)*Dmin
        zt=Hcoil/2-(i)*Dmin
        
        vertices.append([x1,y1,zb])
        vertices.append([x1,y1,zt]) # up the side!
        
        Q1.append([zb,wb[i]])
        Q1.append([zt,wb[i]]) # up the side!     
        
        phi=np.linspace(wb[i],wrev1[i],Nphi)
        # around circumference
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zt])
            Q1.append([Hcoil/2-i*Dmin,phi[j]])
        
        x2=rcoil*np.cos(wrev1[i])
        y2=rcoil*np.sin(wrev1[i])
        zbnext=-Hcoil/2+(i+bt+1)*Dmin
        
        vertices.append([x2,y2,zbnext]) # down the side!
        Q1.append([zbnext,wrev1[i]]) # down the side!
        
        phi=np.linspace(wrev1[i],wb[i+1],Nphi)
        zbnext=-Hcoil/2+(i+bt+1)*Dmin
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zbnext])
            Q1.append([zbnext,phi[j]])
    
    ## Q1-Q2
    offset=np.pi/2
    out1=np.pi/4
    vertices.append([rcoil*np.cos(out1),rcoil*np.sin(out1),zbnext])
    Q1.append([zbnext,out1])
    
    zb=-Hcoil/2+(bt)*Dmin
    B1.append([zbnext,out1])
    B1.append([zb,out1])
    phi=np.linspace(out1,wrev1[0]+offset,Nphi)
    
    for j in range(0,Nphi):
        vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zb])  
        Q2.append([zb,phi[j]])  
        
    #Q2
    for i in range(0,Nt):
        x1=rcoil*np.cos(wrev1[i]+offset) 
        y1=rcoil*np.sin(wrev1[i]+offset)
        zb=-Hcoil/2+(i+bt)*Dmin
        zt=Hcoil/2-(i)*Dmin
        
        vertices.append([x1,y1,zb])
        vertices.append([x1,y1,zt]) # up the side!
        Q2.append([zb,wrev1[i]+offset])
        Q2.append([zt,wrev1[i]+offset]) # up the side!
        
        phi=np.linspace(wrev1[i]+offset,wb[i]+offset,Nphi)
        # around circumference
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zt])
            Q2.append([zt,phi[j]])
        
        x2=rcoil*np.cos(wb[i]+offset)
        y2=rcoil*np.sin(wb[i]+offset)
        zbnext=-Hcoil/2+(i+bt+1)*Dmin
        
        vertices.append([x2,y2,zbnext]) # down the side!
        Q2.append([zbnext,wb[i]+offset]) # down the side!
        
        phi=np.linspace(wb[i]+offset,wb[::-1][i+1]+offset,Nphi)

        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zbnext])
            Q2.append([zbnext,phi[j]])
    
    
    ## Q2-Q3      
    zb=-Hcoil/2+(bt-1)*Dmin
    
    vertices.append([rcoil*np.cos(out1+offset),rcoil*np.sin(out1+offset),zbnext])
    Q2.append([zbnext,out1+offset])
    
    B2.append([zbnext,out1+offset])
    B2.append([zb,out1+offset])
    offset2=np.pi
    phi=np.linspace(out1+offset,wb[0]+offset2,Nphi)
    
    for j in range(0,Nphi):
        vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zb])          
        Q3.append([zb,phi[j]])   
        
        
    #Q3    
    for i in range(0,Nt):
        x1=rcoil*np.cos(wb[i]+offset2) 
        y1=rcoil*np.sin(wb[i]+offset2)
        zb=-Hcoil/2+(i+bt)*Dmin
        zt=Hcoil/2-(i)*Dmin
        
        vertices.append([x1,y1,zb])
        vertices.append([x1,y1,zt]) # up the side!
        Q3.append([zb,wb[i]+offset2])
        Q3.append([zt,wb[i]+offset2]) # up the side!   
     
        phi=np.linspace(wb[i]+offset2,wrev1[i]+offset2,Nphi)
        # around circumference
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zt])
            Q3.append([zt,phi[j]])
        
        x2=rcoil*np.cos(wrev1[i]+offset2)
        y2=rcoil*np.sin(wrev1[i]+offset2)
        zbnext=-Hcoil/2+(i+bt+1)*Dmin
        
        vertices.append([x2,y2,zbnext]) # down the side!
        Q3.append([zbnext,wrev1[i]+offset2]) # down the side!
        
        phi=np.linspace(wrev1[i]+offset2,wb[i+1]+offset2,Nphi)
        
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zbnext])
            Q3.append([zbnext,phi[j]])
            
            

    ## Q3-Q4   
    zb=-Hcoil/2+(bt-1)*Dmin
    vertices.append([rcoil*np.cos(out1+offset2),rcoil*np.sin(out1+offset2),zbnext])
    Q3.append([zbnext,out1+offset2])
    B3.append([zbnext,out1+offset2])
    B3.append([zb,out1+offset2])
    offset3=3*np.pi/2
    phi=np.linspace(out1+offset2,wrev1[0]+offset3,Nphi)
    
    for j in range(0,Nphi):
        vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zb])   
        Q4.append([zb,phi[j]])   
    
    #Q4
    for i in range(0,Nt):
        x1=rcoil*np.cos(wrev1[i]+offset3) 
        y1=rcoil*np.sin(wrev1[i]+offset3)
        zb=-Hcoil/2+(i+bt)*Dmin
        zt=Hcoil/2-(i)*Dmin
        
        vertices.append([x1,y1,zb])
        vertices.append([x1,y1,zt]) # up the side!
        Q4.append([zb,wrev1[i]+offset3])
        Q4.append([zt,wrev1[i]+offset3]) # up the side!
        
        phi=np.linspace(wrev1[i]+offset3,wb[i]+offset3,Nphi)
        # around circumference
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zt])
            Q4.append([zt,phi[j]])
        
        x2=rcoil*np.cos(wb[i]+offset3)
        y2=rcoil*np.sin(wb[i]+offset3)
        zbnext=-Hcoil/2+(i+bt+1)*Dmin
        vertices.append([x2,y2,zbnext]) # down the side!
        Q4.append([zbnext,wb[i]+offset3]) # down the side!
        
        phi=np.linspace(wb[i]+offset3,wb[::-1][i+1]+offset3,Nphi)
        for j in range(0,Nphi):
            vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zbnext])
            Q4.append([zbnext,phi[j]])
    
    ## Q4-Q1   
    zb=-Hcoil/2
    vertices.append([rcoil*np.cos(out1+offset3),rcoil*np.sin(out1+offset3),zbnext])
    Q4.append([zbnext,out1+offset3])
    
    B4.append([zbnext,out1+offset3])
    B4.append([zb,out1+offset3])
    
    phi=np.linspace(out1+offset3,wb[1],Nphi)
    for j in range(0,Nphi):
        vertices.append([rcoil*np.cos(phi[j]),rcoil*np.sin(phi[j]),zb])   
        Q5.append([zb,phi[j]])   
    
    spacing=xwires[1:]-xwires[:-1]
    print('argmin ',np.argmin(spacing))
    print('min wire spacing',np.min(np.abs(spacing)))
    print('min wire spacing',2*xwires[0])
    
    B1=np.array(B1)
    B2=np.array(B2)
    B3=np.array(B3)
    B4=np.array(B4)
    
    Q1=np.array(Q1)
    Q2=np.array(Q2)
    Q3=np.array(Q3)
    Q4=np.array(Q4)
    Q5=np.array(Q5)
    

    Backside=[B1,B2,B3,B4]
    Coils=[Q1,Q2,Q3,Q4,Q5]
    
    cd = source.current.Line(curr=1,vertices=vertices)
    
    # create collection
    c = Collection(cd)
    
    return [c,Coils,Backside]

def MakeQuadDXF(fcoords,rc_,fname):
    
    doc = ezdxf.new('R12')
    doc.units = units.MM
    msp = doc.modelspace()

    doc.layers.new(name='MyLines', dxfattribs={'linetype': 'DASHED', 'color': 7})
    my_lines = doc.layers.get('MyLines')
    
    for i in range(1,len(fcoords)):
        c0=rc_*fcoords[i-1][1]
        z0=fcoords[i-1][0]
        
        ci=rc_*fcoords[i][1]
        zi=fcoords[i][0]
        line = msp.add_line((c0, z0), (ci, zi), dxfattribs={'layer': 'MyLines','lineweight':3})
    #my_lines.dxf.linetype = 'DOTTED'

    doc.saveas(fname+".dxf")
    
        
    return None

def MakeZCoilDXF(fcoords,rc_,fname):
    
    doc = ezdxf.new('R12')
    doc.units = units.MM
    msp = doc.modelspace()

    doc.layers.new(name='MyLines', dxfattribs={'linetype': 'DASHED', 'color': 7})
    my_lines = doc.layers.get('MyLines')
    
    for i in range(1,len(fcoords)):
        c0=rc_*fcoords[i-1][1]
        z0=fcoords[i-1][0]
        
        ci=rc_*fcoords[i][1]
        zi=fcoords[i][0]
        line = msp.add_line((c0, z0), (ci, zi), dxfattribs={'layer': 'MyLines','lineweight':3})
    #my_lines.dxf.linetype = 'DOTTED'

    doc.saveas(fname+".dxf")
    
        
    return None
    

    

