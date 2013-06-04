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

fig3=p.figure()
a3=fig3.add_subplot(1,1,1)
a3.set_xlabel('k [rad/m]')




a3.set_ylabel('KE Tendencies [$10^{'+strexpon+'} m^3/s^3$]')

pic3name='KEslope'


for pref in ['5K','10K','20K']:
    filename=pref+'JET.pck'
    file1=open(filename,'r')
    BDATALOADED=pk.load(file1)
    file1.close()
    [amp0,amp1,amp2,amp3,amp4,amp5,amp6,amp7,amp8,amp9,amp10,ktmp,count,kchoicelist]=BDATALOADED
    istr=0
    print(' amp0 rebuilt = %lf ',sum(amp0[istr:]))
    A=[]

#
# #############################################################################
# #
# # PLOT
# #
#
    lx=30
    ly=15



# figure('units','centimeters','position', ...
#          [0 0 lx ly],'paperpositionmode','auto')
#
#
    dk=ktmp[9]-ktmp[8]
    A +=[cff_scale*amp0[istr:]/dk]
    K=ktmp[istr:]
    LK=K.shape[0]
#number1=K[0]
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


