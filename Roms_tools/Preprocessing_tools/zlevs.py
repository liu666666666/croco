#  
#  This file is copied from ROMSTOOLS
#
#  ROMSTOOLS is free software you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation either version 2 of the License,
#  or (at your option) any later version.
#
#  ROMSTOOLS is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#
#  Copyright (c) 2002-2006 by Pierrick Penven 
#  e-mail:Pierrick.Penven@ird.fr  
#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import numpy as np
from math import cosh, exp,sinh, tanh

def csf(sc, theta_s,theta_b):
    if (theta_s > 0 ):
        csrf=(1-cosh(sc*theta_s))/(cosh(theta_s)-1)
    else:
        csrf=-sc**2
    if (theta_b > 0):
        h = (exp(theta_b*csrf)-1)/(1-exp(-theta_b))
    else:
        h  = csrf
    return h



################################################################

def zlevs(h,zeta,theta_s,theta_b,hc,N,typestring,vtransform=1):
    '''
    function z = zlevs(h,zeta,theta_s,theta_b,hc,N,typestring,vtransform)

  function z = zlevs(h,zeta,theta_s,theta_b,hc,N,vtransform,typestring)

  this function compute the depth of rho or w points for ROMS

  On Input:

    typestring    'r': rho point 'w': w point 
    vtransform  1=> old v transform  2=>new v transform 
#  On Output:
#
#    z       Depths (m) of RHO- or W-points (3D matrix).
# 
#  Further Information:  
#  http://www.brest.ird.fr/Roms_tools/

#  
#  This file is copied from ROMSTOOLS
#
#  ROMSTOOLS is free software you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation either version 2 of the License,
#  or (at your option) any later version.
#
#  ROMSTOOLS is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#
#  Copyright (c) 2002-2006 by Pierrick Penven 
#  e-mail:Pierrick.Penven@ird.fr  
#
'''
    
    if vtransform==1:
        print(['WARNING no vtransform defined']) 
        print(['Default S-coordinate system use : Vtransform=1 (old one)'])
    
    #
    # Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
    
    
    sc_r=np.zeros(N)
    Cs_r=np.zeros(N)
    sc_w=np.zeros(N+1)
    Cs_w=np.zeros(N+1)
    temp=np.arange(-N+1,1,1)
    if (vtransform == 2):
        ds=1./N
        if typestring=='w':
    
          sc_w[0] = -1.0
          sc_w[N] =  0
          Cs_w[0] = -1.0
          Cs_w[N] =  0
          
          sc_w[1:N-1] = ds*temp
          Cs_w=csf(sc_w, theta_s,theta_b)
          N=N+1
    
    #    print(['===================================='])
    #    for k=N:-1:1
    #        print(['Niveau S=',num2str(k),' Cs=',num2str( Cs_w(k), '#8.7f')])
    #    
    #    print(['===================================='])
    
        else:
    
          sc= ds*(temp-0.5)    
          Cs_r=csf(sc, theta_s,theta_b)
          sc_r=sc
    #    print(['===================================='])
    #    for k=N:-1:1
    #        print(['Niveau S=',num2str(k),' Cs=',num2str( Cs_r(k), '#8.7f')])
    #    
    #    print(['===================================='])
        
    
    else:
        cff1=1./sinh(theta_s)
        cff2=0.5/tanh(0.5*theta_s)
        if typestring=='w':
            sc=np.arange(-N,1,1)/float(N)
            N=N+1
        else:
            sc=np.arange(-N+0.5,0,1)/float(N)
        
        Cs=(1.-theta_b)*cff1*np.sinh(theta_s*sc)+theta_b*(cff2*np.tanh(theta_s*(sc+0.5))-0.5)
    #    print(['===================================='])
    #    for k=N:-1:1
    #        print(['Niveau S=',num2str(k),' Cs=',num2str( Cs(k), '#8.7f')])
    #    
    #    print(['===================================='])
    
    #
    # Create S-coordinate system: based on model topography h(i,j),
    # fast-time-averaged free-surface field and vertical coordinate
    # transformation metrics compute evolving depths of of the three-
    # dimensional model grid.
    #  
    hinv=1/h
    if type(h) is not np.ndarray:
        M=L=1
    else:
        (M,L)=h.shape
    z=np.zeros((N,M,L))
    if (vtransform == 2):
        if typestring=='w':
            cff1=Cs_w
            cff2=sc_w+1
            sc=sc_w
        else:
            cff1=Cs_r
            cff2=sc_r+1
            sc=sc_r
        
        h2=(h+hc)
        cff=hc*sc
        h2inv=1./h2
        for k in range(N):
            z0=cff[k]+cff1[k]*h
            z[k,:,:]=z0*h/(h2) + zeta*(1+z0*h2inv)
        
    else:
        cff1=Cs
        cff2=sc+1
        cff=hc*(sc-Cs)
        cff2=sc+1
        for k in range(N):
            z0=cff[k]+cff1[k]*h
            z[k,:,:]=z0+zeta*(1.+z0*hinv)
    
    
    return z
            
        

