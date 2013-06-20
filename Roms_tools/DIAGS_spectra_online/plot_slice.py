# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 11:30:18 2013

@author: soufflet

program to plot 3D slices and create an animation.

"""

import sys
from pylab import *
import pyroms
from utilities import cut_var

def main(argv):
    print("Hola!!!")
    filename=argv[0]
    level=int(argv[2])
    varname=argv[1]
    print('reading netcdf file: %s',filename)  
    lims=[1,100,1,400]
    ncfile=pyroms.io.Dataset(filename)
    var1=ncfile.variables[varname][-100:,level,:,:]
    var1=cut_var(var1.transpose(),lims)
    X=ncfile.variables['x_rho'][:,:]
    Y=ncfile.variables['y_rho'][:,:]
    X=cut_var(X.transpose(),lims)
    Y=cut_var(Y.transpose(),lims)
        
    for n in range(100):
        newfilename='slices/'+varname+'_'+str(level)+'.{num:04d}.png'.format(num=n)
        #from mpl_toolkits.mplot3d import Axes3D
        fig=figure()
        ax = fig.add_subplot(111)
        #poly = matplotlib.collections.PolyCollection([var,var])
        #ax.add_collection3d(poly, zs=range(2), zdir='z')
        #im=ax.pcolor(X/1000,Y/1000,var[:,:,0])
        #fig.colorbar(im)
        essai=ax.pcolor(X/1000,Y/1000,var1[:,:,n])
        essai.set_clim([var1[:,:,0].min(),var1[:,:,0].max()])    
        fig.colorbar(essai)
        savefig(newfilename)
    #ax.draw()
        #axis([-3,3,-3,3])       
#nshow()
    #
    #show()
    #plt.show()
#def plot_withwidget(var1,time,X,Y,depth):
#    from matplotlib.widgets import Slider, Button
#    fig1=figure()
#    subplots_adjust(left=0.25, bottom=0.25)
#
#    ax=fig1.add_subplot(111)    
#    axcolor = 'lightgoldenrodyellow'
#    axtime = fig1.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
#    resetax = fig1.add_axes([0.8, 0.025, 0.1, 0.04])
#    
#    
#    
#    stime = Slider(axtime, 'time',0, time.shape[0],valfmt='%d', valinit=0)
#    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
#    essai=ax.pcolor(X/1000,Y/1000,var1[:,:,0])
#    essai.set_clim([var1[:,:,0].min(),var1[:,:,0].max()])    
#    fig1.colorbar(essai)
#    fig1.suptitle('Vertical velocity ')   
#    
##    axamp  = axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
#    
#   
##    samp = Slider(axamp, 'Amp', 0.1, 10.0, valinit=a0)
##    ax1=fig1.axes[2]
##    labels=[line.get_label() for line in ax1.lines]
##    fontP = matplotlib.font_manager.FontProperties()
##    fontP.set_size('small')
##    fig1.legend(ax1.lines,labels,'center left',prop=fontP)
##    fig1.canvas.widgetlock(axfreq)
#    
#    def update(val):
#        freq = int(sfreq.val)
##        Dsub=[d for d in D if d[1]==freq]
#        #fig1.clear()
#        plot_profile(fig1,S[freq,1:],T[freq,1:],mini,maxi,depthlist)
#        pylab.draw()
#        
#        
#    sfreq.on_changed(update)
#    
#   
#    
##    fig1.canvas.widgetlock(resetax)
##    resetax.figure.canvas.widgetlock(resetax)
#    button.ax.figure.canvas.widgetlock(button)
##    sfreq.ax.figure.canvas.widgetlock(sfreq)
#    def reset(event):
#        sfreq.reset()
##        just_plot(fig1,[f0],Dsub,suf,pathname)
##        draw()
#    button.on_clicked(reset)

if __name__ == "__main__":
    #print("Hello ..."+sys.argv[1:])
    main(sys.argv[1:])
    bashCommand = "ls slices  "
    import os
    os.system(bashCommand)

