import numpy as np
import pylab as p
import pickle as pk

cff_scale=1.e7
number1=cff_scale
from math import log10 
expon=-int(np.floor(log10(abs(number1))))
coef=number1*10.**expon
print expon,coef
#cff_scale=10.**4
strexpon=str(expon)

for pref in ['5K','10K','20K']:
    resol=pref
    filename=pref+'JET.pck'
    file1=open(filename,'r')
    BDATALOADED=pk.load(file1)
    file1.close()
    [amp0,amp1,amp2,amp3,amp4,amp5,amp6,amp7,amp8,amp9,amp10,ktmp,count,kchoicelist]=BDATALOADED
    istr=0
    print(' amp0 rebuilt = %lf ',sum(amp0[istr:]))
    A=[]

    # # PLOT
    # #
    #
    lx=30
    ly=15

    dk=ktmp[9]-ktmp[8]

    K=ktmp[istr:]
    LK=K.shape[0]

    A +=[cff_scale*amp0[istr:]/dk]
    A+=[cff_scale*amp2[istr:]/dk]
    A+=[cff_scale*amp3[istr:]/dk]
    A+=[cff_scale*amp4[istr:]/dk]
    A+=[cff_scale*amp5[istr:]/dk]
    A+=[cff_scale*amp6[istr:]/dk]
    A+=[cff_scale*amp7[istr:]/dk]
    A+=[cff_scale*amp8[istr:]/dk]
    A+=[cff_scale*amp9[istr:]/dk]
    A+=[cff_scale*amp10[istr:]/dk]
    amp11=amp4-(amp2+amp3+amp5+amp6+amp7+amp10)
    A+=[cff_scale*amp11[istr:]/dk]

    fig=p.figure()
    a=fig.add_subplot(1,1,1)
    a.plot(K,A[2],'b',label='horizontal diffusion')
    a.plot(K,A[1],'g',label='horizontal advection')
    a.plot(K,A[3],'r',label=r'$\frac{\partial{u}}{\partial{t}}$')
    a.plot(K,A[4],'c',label='P')
    a.plot(K,A[5],'b--',label='vertical diffusion',)
    a.plot(K,A[6],'g--',label='vertical advection')
    a.plot(K,A[9],'r--',label='Coriolis')
    a.set_xscale('log')
    a.legend()
    pic1name=resol+'budget'
    fig.suptitle(resol+' spectral budget')
    fig.savefig(pic1name)

    fig2=p.figure()
    a2=fig2.add_subplot(1,1,1)
    a2.plot(K,A[7],'c--',label='Residual')
    a2.plot(K,A[10],'k-*',label='Res(spectral)')
    a2.set_xscale('log')
    a2.legend()
    strexpon=str(expon)
    a2.set_xlabel('k [rad/m]')
    a2.set_ylabel('KE Tendencies [$10^{'+strexpon+'} m^3/s^3$]')
    pic2name=resol+'residual'
    fig2.suptitle(resol+' residual_pt')
    fig2.savefig(pic2name)

    fig3=p.figure()
    a3=fig3.add_subplot(1,1,1)
    a3.set_xlabel('k [rad/m]')
    a3.set_ylabel('KE Tendencies [$10^{'+strexpon+'} m^3/s^3$]')
    pic3name=resol+'KEslope'


    a3.plot(K,amp0/dk,label=pref)

    a3.plot(K,K**(-2)/cff_scale,'r',label=r'$k^{-2}$')
    a3.plot(K,K**(-5./3)/cff_scale,'b--',label=r'$k^{-5/3}$')
    a3.plot(K,K**(-3)/cff_scale,'c+',label=r'$k^{-3}$')
    a3.set_yscale('log')
    a3.set_xscale('log')
    a3.set_xlabel('Wavenumber k [rad/m]')
    a3.set_title('KE spectrum')
    a3.legend()
    fig3.savefig(pic3name)
p.show()


