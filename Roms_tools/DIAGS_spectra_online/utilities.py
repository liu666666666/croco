# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 16:53:58 2013

@author: soufflet
"""

import numpy as np
class gridinfo:
    '''
    Class redifning the gridinfo structure
    '''
    id      = '';
    name    = 'baroclinic jet grid 5km'
    grdfile = 'test.nc'
    N       = 30
    thetas  = 5
    thetab  = 0
    hc      = 400
    Method  = 1

    def setId(self,idname):
        self.id=idname



def cut_var(var,lims):
    '''
    function to extract a subset of a variable
    '''
    minL=lims[0]
    maxL=lims[1]
    minM=lims[2]
    maxM=lims[3]
    if var.ndim==3:
              var=np.squeeze(var[minL:maxL,minM:maxM,:])
    elif var.ndim==2:
              var=np.squeeze(var[minL:maxL,minM:maxM])
    else:
              print('wrong var size in CUT_VAR: #d',var.ndim)
    return var

def partialx(var,pm):
    '''
    compute the relative x derivative of var, pm is the metrics
    array of var.size [Mp,Lp,1]. This is done for a 2d or 3d variable
    '''
    (Mp,Lp)=pm.shape
    if var.ndim==2:
        (Mv,Lv)=var.shape
        dxvar=np.zeros((Mv,Lv))
        dxvar[1:-1,:]=0.5*pm[1:-1,:]*(var[2:,:]-var[:-2,:])
        dxvar[0,:]=dxvar[1,:]
        dxvar[-1,:]=dxvar[-2,:]
    elif var.ndim==3:  # 3d variable
        (Mv,Lv,Kv)=var.shape
        dxvar=np.zeros((Mv,Lv,Kv))
        pm2=pm.reshape((Mp,Lp,1))
        dxvar[1:-1,:,:]=0.5*pm2[1:-1,:,0]*(var[2:,:,:]-var[:-2,:,:])
        dxvar[0,:,:]=dxvar[1,:,:]
        dxvar[-1,:,:]=dxvar[-2,:,:]
    elif var.ndim==4: # time variable
        (Mv,Lv,Kv,Tv)=var.shape
        dxvar=np.zeros((Mv,Lv,Kv,Tv))
        pm2=pm.reshape((Mp,Lp,1,1))
        dxvar[1:-1,:,:,:]=0.5*pm[1:-1,:,0,0]*(var[2:,:,:,:]-var[:-2,:,:,:])
        dxvar[0,:,:,:]=dxvar[1,:,:,:]
        dxvar[-1,:,:,:]=dxvar[-2,:,:,:]
    return dxvar

def partialy(var,pn):
    '''
    compute the relative y derivative of var, pn is the metrics
    array of var.size [Mp,Lp,1]. This is done for a 2d or 3d variable
    '''
    varx=var.swapaxes(0,1)
    pnx=pn.swapaxes(0,1)
    dxvar=partialx(varx,pnx)
    dyvar=dxvar.swapaxes(0,1)
    return dyvar

def rho2d(var,naxis=1):
    u=var.swapaxes(0,naxis)
    (M,P)=u.shape
    u2rho=np.zeros((M+1,P))
    u2rho[1:-1,:]=0.5*(u[:-1,:]+u[1:,:])
    u2rho[0,:]=u2rho[1,:]
    u2rho[-1,:]=u2rho[-2,:]
    return u2rho.swapaxes(0,naxis).transpose()



# These window functions are similar to those found in the Windows toolbox of MATLAB
# Note that numpy has a couple of Window functions already:
# See: hamming, bartlett, blackman, hanning, kaiser
#copied from http://leohart.wordpress.com/2006/01/29/hello-world/

def tukeywin(window_length, alpha=0.5):
    '''The Tukey window, also known as the tapered cosine window, can be regarded as a cosine lobe of width \alpha * N / 2
    that is convolved with a rectangle window of width (1 - \alpha / 2). At \alpha = 1 it becomes rectangular, and
    at \alpha = 0 it becomes a Hann window.

    We use the same reference as MATLAB to provide the same results in case users compare a MATLAB output to this function
    output

    Reference
    ---------
    http://www.mathworks.com/access/helpdesk/help/toolbox/signal/tukeywin.html

    '''
    # Special cases
    if alpha <= 0:
        return np.ones(window_length) #rectangular window
    elif alpha >= 1:
        return np.hanning(window_length)

    # Normal case
    x = np.linspace(0, 1, window_length)
    w = np.ones(x.shape)

    # first condition 0 <= x < alpha/2
    first_condition = x<alpha/2
    w[first_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[first_condition] - alpha/2) ))

    # second condition already taken care of

    # third condition 1 - alpha / 2 <= x <= 1
    third_condition = x>=(1 - alpha/2)
    w[third_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[third_condition] - 1 + alpha/2)))

    return w

def integ_fft2d(tmpamp,dx):
    (L,M)=tmpamp.shape    
    Lh=L/2.
    Mh=M/2.
    cff=2*np.pi/dx
    i2L=1./L**2
    i2M=1./M**2
    a=np.zeros((L,1))
    a[:,0]=np.arange(L)
    
    kk=np.tile(a-Lh,(1,M))
    
    b=np.zeros((1,M))
    b[0,:]=np.arange(M)
    ll=np.tile(b-Mh,(L,1))
    Kh=i2L*kk**2+i2M*ll**2
    Kh=cff*np.sqrt(Kh)
    #Kh=Kh.transpose();
    #
    # Ring integration method with non constant dk
    # accounts for different lengths in xi and eta
    #
    if M>L:
        dcase=1
        Lmin=L
        Lmax=M  #xi dim is bigger
    else:
        dcase=0
        Lmin=M
        Lmax=L
    leng=Lmin/2
    #+int(Lmax/Lmin)-1
    amp=np.zeros(leng)
    ktmp=np.zeros(leng)
    count=np.zeros(leng)
    dk=np.zeros(leng)
      ### first sample the long wavelengths permitted in longest 
      ### direction with wavenumber 0 in shortest direction
    indaux=0
    for ind in range(int(Lmax/Lmin)-1):
        amp[ind]=tmpamp[Lh+(1-dcase)*(ind+1),Mh+dcase*(ind+1)]
        count[ind]=1
        ktmp[ind]=cff/Lmax*(ind+1)
        dk[ind]=cff/Lmax
        indaux=indaux+1;
      
      ### second do the regular summation over rings with 
      ### step that is function of k step for shortest dimension 
      ### ie we use the largest step.
    
    for ind in range(int(Lmax/Lmin)-1,int(Lmin/2)):
        ktmp[ind]=cff/Lmin*(ind+1-indaux)
    #        print('ktmp: %lf',ktmp[ind])
        II=tmpamp[abs(Kh-ktmp[ind])<cff/(2*Lmin)]
    #        print('II: %d',len(II))
        #II=find(abs(Kh-ktmp(ind))<cff/(2*Lmin))
        #[xI,yI]=find(abs(Kh-ktmp(ind))<cff/(2*Lmin))
    #        print tmpamp[II]
        amp[ind]=II.sum()
    #        print('anp: %lf',amp[ind])
        count[ind]=II.shape[0]
        dk[ind]=cff/Lmin

  
#  figure;imagesc(zob);shading flat;title(['should be nan']);
      ### third correct for spectral defects due to summation over 
      ### finite number of points that imperfectly sample the wavenumber 
      ### space.
    
    cff=2*dk[:]*ktmp[:]/(ktmp[-1]**2)*sum(count)/count[:]
    amp=amp*cff;
      #amp=modif_amp(amp,count,ktmp,dk);
    
    return amp,count,ktmp,dk