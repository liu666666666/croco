import numpy as np
import pickle as pk
import pylab as p
    
import pycomodo as pc

from utilities import rho2d_xperiodic1,extend_torhoshape, tukeywin,integ_fft2d, get_zw_fromnc, horizontal_advection_C4


def compute_budget(resol,depth_integration,tlength):
    ''' Setting main parameters of for the calculation:
        * if mean =1 , we remove the mean of each variable and use the fluctuation instead u'=u-ubar 
        * window =1 we use windowing (useful for non-periodic solution)
        * cff_tukey=1, if windowing is used, uses the Tukey method
        * salinity=1 if we uses salinity in density calculation
        * wind=1 if wind forcing is used'''

    mmean=0
    window=0
    cff_tukey=1
    salinity=0
    wind=0


    ''' Path of the output files, resolution, box set, vertical level list, iterations... '''
    #resol='5K'
    nres=resol.split('K')[0]
    dx=float(nres)*1000
    Nx=int(500000/dx)
    Ny=int(2000000/dx)

    lims=[0,Nx,0,Ny]
    #root='/data/models/JET/Last3Months/'+resol #+'JET/'
    root='/data/models/JET/LastYear/'+resol#+'JET/AsselinFrom18M/'#+resol
    model='jet_Lastyear'
    #model='jet'#_100Days'

    '''physical constants and plotting parameters'''

    g=9.81
    rho0=1024
    cff_scale=1.e7
    xmin=0.9e-5
    xmax=1.e-3
    ymin=-1
    ymax=1  # limits for plots

    historyfilename=root+model+'_avg.nc'
    print historyfilename

    print('reading netcdf file: %s',historyfilename)
    nchis=pc.Archive(historyfilename)
        
    #overload grid_info function
    zw=get_zw_fromnc(nchis).squeeze()
    N=zw.shape[0]-1
    zlayers=0.5*(zw[:-1]+zw[1:])
    level=np.argmin(np.abs(zlayers+depth_integration))
    kchoicelist=range(N)[level:]
    kchoice=level
    print('Integrating from depth: ',zlayers[kchoice])
    dz=np.zeros(zlayers.shape)
    dz=zw[1:]-zw[:-1];dz=dz.squeeze()
    dz_klist=dz[kchoicelist]
    dz=dz[kchoicelist]
    depthmax=int(zlayers[kchoice])
    
    timevar=nchis.variables['time']
    totaliterations=timevar.shape[0]
    #it2=totaliterations-tlength-1
    #it1=totaliterations-tlength-120
    it2=totaliterations-1
    it1=totaliterations-tlength

    suffixe='_averaged_at'+str(depthmax)+'m_fromday_'+str(tlength)
    
    #loaded=True
    #save=False
    #
    loaded=False
    save=True

    filename=resol+'JET'+suffixe+'.pck'


    ''' Different terms are already computed, online=1, or needs to be calculated'''
    online_T=1
    online_P=1
    online_Cor=1
    online_Dv=1

    nlevels=len(kchoicelist)

    
    
    
    tlength=it2-it1+1
    print ("iterations:",it1,it2,tlength)
    if not loaded:
        
        diagfilename =root+model+'_diaM_avg.nc'
        srhoshape=nchis.variables['s_rho'].shape
        print diagfilename
        print("srho shape:",srhoshape)
        print('reading netcdf file: %s',diagfilename) 
        ncdiags=pc.Archive(diagfilename)

        print('domain limits: ',lims) 

        U=nchis.variables['u']
        V=nchis.variables['v']
        pm=nchis.variables['pm']
        pn=nchis.variables['pn']
        (Mu,Lu)=U.shape[-2:]
        (Mv,Lv)=V.shape[-2:]
        M=Mu
        L=Lv

        # ###############################################################
        #
        print("initialistion of sum variables:")
        e_rate=0
        e_hdis=0
        e_hadv=0
        e_hdpr=0
        e_zadv=0
        e_zdis=0
        e_cor=0
        e_res=0
        e_wb=0

        print("initialisation des var de calcul de spectre")

        tmpamp0=np.zeros((M,L))
        tmpamp1=np.zeros((M,L))
        tmpamp2=np.zeros((M,L))
        tmpamp3=np.zeros((M,L))
        tmpamp4=np.zeros((M,L))
        tmpamp5=np.zeros((M,L))
        tmpamp6=np.zeros((M,L))
        tmpamp7=np.zeros((M,L))
        tmpamp8=np.zeros((M,L))
        tmpamp9=np.zeros((M,L))
        tmpamp10=np.zeros((M,L))
        tmpamp22=np.zeros((M,L))

        if window==1:
            cff_tukey=0.25
            wdw1=tukeywin(M,-1.0)
            #wdw1=np.ones(L)
            wdw2=tukeywin(L,cff_tukey)
            #print(wdw1.shape)
            #print(wdw2.shape)
            wdw=np.zeros((wdw1.shape[0],wdw2.shape[0]))
            #wdw=np.dot(wdw1.transpose(),wdw2)
            wdw=wdw1.reshape(wdw1.shape[0],1)*wdw2.reshape(1,wdw2.shape[0])


##
##        C4_hadvection_u=np.zeros((tlength,nlevels,Mu,Lu))
##        C4_hadvection_v=np.zeros((tlength,nlevels,Mv,Lv))
##
##        UP3_hadvection_u=ncdiags.variables['u_xadv'][-tlength:,kchoice:,:,:]
##        UP3_hadvection_u=UP3_hadvection_u+ncdiags.variables['u_yadv'][-tlength:,kchoice:,:,:]
##
##        UP3_hadvection_v=ncdiags.variables['v_xadv'][-tlength:,kchoice:,:,:]
##        UP3_hadvection_v=UP3_hadvection_v+ncdiags.variables['v_yadv'][-tlength:,kchoice:,:,:]
##

        ik=-1

        #
        # VERTICAL LOOP ================================
        #

        for kchoice in kchoicelist:
            ik=ik+1
            
            #
            # TEMPORAL LOOP ================================
            #

            for it in np.arange(it1,it2):
                '''
                In ROMS history files contain initialisation fields (t=0) 
                so start to read from iteration number 1
                '''
                ithis=it#+1

                print('Treating rec #d',it)
                timehis=nchis.variables['scrum_time'][ithis]
                timediaM=ncdiags.variables['scrum_time'][it]
                print('which, in his file, correspond to: ', int(timehis))
                print('which, in diaM file, correspond to: ', int(timediaM))

                levelhis=nchis.variables['s_rho'][kchoice]
                leveldiaM=ncdiags.variables['s_rho'][kchoice]
                print('which, in his file, correspond to :', levelhis)
                print('Depth: ',zw[kchoice])
                print('which, in diaM file, correspond to :', leveldiaM)

                u=nchis.variables['u'][ithis,kchoice,:,:]
                v=nchis.variables['v'][ithis,kchoice,:,:]
                Hz=dz[ik]
#                advu1,advv1=horizontal_advection_C4(u,v,Hz,pm,pn)

#                C4_hadvection_u[it,ik,:,:]=advu1
#                C4_hadvection_v[it,ik,:,:]=advv1

                u=rho2d_xperiodic1(u)
                v=rho2d_xperiodic1(v,naxis=0)

        #### UPSTREAM BIASED NL TERMS ####
        #### computed as if in the rhs (- sign already included).

        #==============================================================================
                advu2=ncdiags.variables['u_xadv'][it,kchoice,:,:]
                advu2=advu2+ncdiags.variables['u_yadv'][it,kchoice,:,:]
                

                advv2=ncdiags.variables['v_xadv'][it,kchoice,:,:]
                advv2=advv2+ncdiags.variables['v_yadv'][it,kchoice,:,:]

                advu2[:,-1]= advu2[:,0]
                advv2[:,0]= advv2[:,-2]
                advv2[:,-1]= advv2[:,1]

                advu2=rho2d_xperiodic1(advu2)
                advv2=rho2d_xperiodic1(advv2,naxis=0)
 
                advu1=advu2
                advv1=advv2

        #
        # #### implicit dissipation
        # #### computed as in the rhs
                dissu=advu2-advu1
                dissv=advv2-advv1
        #
        # #### explicit dissipation
        # #### computed as if in the rhs (- sign already included).

    ## dissu2=ncdiags.variables['u_hmix'][it,kchoice,:,:]
    ## dissu2=rho2d_xperiodic1(dissu2)
    ## dissu2=cut_var(dissu2,lims)
    ##
    ## dissv2=ncdiags.variables['v_hmix'][it,kchoice,:,:]
    ## dissv2=rho2d_xperiodic1(dissv2,naxis=0)
    ## dissv2=cut_var(dissv2,lims)
    ##
    ## dissu=dissu2
    ## dissv=dissv2
    ##
        # #==============================================================================
                if wind:
                    print("let's not do WIND right now!!!")
                    exit

        #
        # #### temporal derivative
                if online_T==1:
                    dudt=ncdiags.variables['u_rate'][it,kchoice,:,:]
                    dvdt=ncdiags.variables['v_rate'][it,kchoice,:,:]

                    dudt[:,-1]= dudt[:,0]
                    dvdt[:,0]= dvdt[:,-2]
                    dvdt[:,-1]= dvdt[:,1]

                    dudt=rho2d_xperiodic1(dudt)
                    dvdt=rho2d_xperiodic1(dvdt,naxis=0)
                else:
                    print("let's not do offline_T right now!!!")
                    exit


        #
        # #### wb term
                cffb=-g/rho0 # conversion from density to buoyancy
                w=nchis.variables['w'][ithis,kchoice,:,:]
                w[:,-1]=w[:,1]
 
                w=w-w.mean()

                if salinity:
                    print("let's not do salinity right now!!!")
                    exit
                    T=nchis.variables['temp'][ithis,kchoice,:,:]
                    S=nchis.variables['salt'][ithis,kchoice,:,:]
                    #rho=rho_potential(T,S)
                else:
                    rho=nchis.variables['temp'][ithis,kchoice,:,:]

                b=cffb*(rho-rho.mean()) 

        #
        # #### pressure term
                if online_P==1:
                    dxp=ncdiags.variables['u_Prsgrd'][it,kchoice,:,:]
                    dyp=ncdiags.variables['v_Prsgrd'][it,kchoice,:,:]

                    dxp[:,-1]= dxp[:,0]
                    dyp[:,0]=dyp[:,-2]
                    dyp[:,-1]=dyp[:,1]
                    dxp=rho2d_xperiodic1(dxp)
                    dyp=rho2d_xperiodic1(dyp,naxis=0)
                else:
                    print("let's not do that right now!!!")
                    exit
        #
        # # Coriolis
        # #
                if online_Cor==1:
                    ucor=ncdiags.variables['u_cor'][it,kchoice,:,:]
                    vcor=ncdiags.variables['v_cor'][it,kchoice,:,:]
                    #fix diag var
                    ucor[:,-1]= ucor[:,0]
                    vcor[:,0]= vcor[:,-2]
                    vcor[:,-1]= vcor[:,1]

                    ucor=rho2d_xperiodic1(ucor)
                    vcor=rho2d_xperiodic1(vcor,naxis=0)
                else:
                    print("let's not do offline_Coriolis right now!!!")
                    exit
        #
        # #### vertical dissipation
                if online_Dv==1:
                    ukpp=ncdiags.variables['u_vmix'][it,kchoice,:,:]
                    vkpp=ncdiags.variables['v_vmix'][it,kchoice,:,:]
                    ukpp[:,-1]= ukpp[:,0]
                    vkpp[:,0]= vkpp[:,-2]
                    vkpp[:,-1]= vkpp[:,1]
                   
                    ukpp=rho2d_xperiodic1(ukpp)
                    vkpp=rho2d_xperiodic1(vkpp,naxis=0)

                else:
                    print("let's not do that right now!!!")
                    exit

        #
        # #### vertical advection
                zadvu=ncdiags.variables['u_vadv'][it,kchoice,:,:]
                zadvv=ncdiags.variables['v_vadv'][it,kchoice,:,:]
                zadvu[:,-1]= zadvu[:,0]
                zadvv[:,0]= zadvv[:,-2]
                zadvv[:,-1]= zadvv[:,1]

                zadvu=rho2d_xperiodic1(zadvu)
                zadvv=rho2d_xperiodic1(zadvv,naxis=0)
        #
        # #### residual
                resu=-ncdiags.variables['u_rate'][it,kchoice,:,:]+\
                    ncdiags.variables['u_xadv'][it,kchoice,:,:]+\
                    ncdiags.variables['u_yadv'][it,kchoice,:,:]+\
                    ncdiags.variables['u_cor'][it,kchoice,:,:]+\
                    ncdiags.variables['u_vmix'][it,kchoice,:,:]+\
                    ncdiags.variables['u_hmix'][it,kchoice,:,:]+\
                    ncdiags.variables['u_vadv'][it,kchoice,:,:]+\
                    ncdiags.variables['u_Prsgrd'][it,kchoice,:,:]
                resu=-resu
                resu[:,-1]= resu[:,0]

                resv=-ncdiags.variables['v_rate'][it,kchoice,:,:]+\
                    ncdiags.variables['v_xadv'][it,kchoice,:,:]+\
                    ncdiags.variables['v_yadv'][it,kchoice,:,:]+\
                    ncdiags.variables['v_cor'][it,kchoice,:,:]+\
                    ncdiags.variables['v_vmix'][it,kchoice,:,:]+\
                    ncdiags.variables['v_hmix'][it,kchoice,:,:]+\
                    ncdiags.variables['v_vadv'][it,kchoice,:,:]+\
                    ncdiags.variables['v_Prsgrd'][it,kchoice,:,:]

                resv=-resv
                resv[:,0]= resv[:,-2]
                resv[:,-1]= resv[:,1]

                resu=rho2d_xperiodic1(resu)
                resv=rho2d_xperiodic1(resv,naxis=0)
          
                if wind and kchoice==kchoicelist[0]:
                    print("let's not do WIND right now!!!")
                    exit

                if mmean==1:
                    # b and w already done
                    u=u-u.mean()
                    v=v-v.mean()
                    advu1=advu1-mean(mean(advu1))
                    advv1=advv1-mean(mean(advv1))
                    advu2=advu2-advu2.mean()
                    advv2=advv2-advv2.mean()
                    dissu=dissu-dissu.mean()
                    dissv=dissv-dissv.mean()
                    dudt=dudt-dudt.mean()
                    dvdt=dvdt-dvdt.mean()
                    dxp=dxp-dxp.mean()
                    dyp=dyp-dyp.mean()
                    ukpp=ukpp-ukpp.mean()
                    vkpp=vkpp-vkpp.mean()
                    zadvu=zadvu-zadvu.mean()
                    zadvv=zadvv-zadvv.mean()
                    ucor=ucor-ucor.mean()
                    vcor=vcor-vcor.mean()
                    resu=resu-resu.mean()
                    resv=resv-resv.mean()
        #
                if window==1:
                    u=u*wdw
                    v=v*wdw
                    #advu1=advu1*wdw
                    #advv1=advv1*wdw
                    advu2=advu2*wdw
                    advv2=advv2*wdw
                    dissu=dissu*wdw
                    dissv=dissv*wdw
                    dudt=dudt*wdw
                    dvdt=dvdt*wdw
                    dxp=dxp*wdw
                    dyp=dyp*wdw
                    ukpp=ukpp*wdw
                    vkpp=vkpp*wdw
                    zadvu=zadvu*wdw
                    zadvv=zadvv*wdw
                    ucor=ucor*wdw
                    vcor=vcor*wdw
                    resu=resu*wdw
                    resv=resv*wdw
                    b=b*wdw
                    w=w*wdw


                #sum all terms to check mean value before projection
                e_rate+=np.mean(u*dudt+v*dvdt)/tlength
                e_hdis+=np.mean(u*dissu+v*dissv)/tlength
                e_hadv+=np.mean(u*advu2+v*advv2)/tlength
                e_hdpr+=np.mean(u*dxp  +v*dyp)/tlength
                e_zadv+=np.mean(u*zadvu+v*zadvv)/tlength
                e_zdis+=np.mean(u*ukpp+v*vkpp)/tlength
                e_cor+=np.mean(u*ucor+v*vcor)/tlength
                e_wb+=np.mean(w*b)/tlength
                e_res+=(e_hadv+e_hdpr+e_zadv+e_zdis+e_cor)/tlength

        #
        #  SPECTRAL PROJECTION
        #
        


        #
                #print('shape 5 (%d,%d)',tmpamp0.max())
                cff=1./((L*M)**2)
                cff1=cff*dz_klist[ik]/sum(dz_klist)/tlength
                fcoefu=np.fft.fft2(u)
                fcoefv=np.fft.fft2(v)
                tmpamp0+=cff1*np.real(np.conj(fcoefu)*fcoefu+np.conj(fcoefv)*fcoefv)        #  0  kinetic energy
        #
                fcoeftmpu=np.fft.fft2(advu1)
                fcoeftmpv=np.fft.fft2(advv1)
                tmpamp1+=cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  1  C4 horiz advection

                fcoeftmpu=np.fft.fft2(advu2)
                fcoeftmpv=np.fft.fft2(advv2)
                tmpamp2+=cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  2  UP3 horiz advection

                fcoeftmpu=np.fft.fft2(dissu)
                fcoeftmpv=np.fft.fft2(dissv)
                tmpamp3+=cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  3  horiz dissipation

                fcoeftmpu=np.fft.fft2(dudt)
                fcoeftmpv=np.fft.fft2(dvdt)
                tmpamp4+=cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  4  time rate of change

                fcoeftmpu=np.fft.fft2(dxp)
                fcoeftmpv=np.fft.fft2(dyp)
                tmpamp5+=cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  5  pressure gradient

                fcoeftmpu=np.fft.fft2(ukpp)
                fcoeftmpv=np.fft.fft2(vkpp)
                tmpamp6+=cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  6  vert dissipation

                fcoeftmpu=np.fft.fft2(zadvu)
                fcoeftmpv=np.fft.fft2(zadvv)
                tmpamp7+=cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  7  vert advection

                fcoeftmpu=np.fft.fft2(resu)
                fcoeftmpv=np.fft.fft2(resv)
                tmpamp8+=cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  8  residual

                fcoeftmpu=np.fft.fft2(w)
                fcoeftmpv=np.fft.fft2(b)
                tmpamp9+=cff1*np.real(np.conj(fcoeftmpu)*fcoeftmpv)                        #  9  baroclinic conversion wb

                fcoeftmpu=np.fft.fft2(ucor)
                fcoeftmpv=np.fft.fft2(vcor)
                tmpamp10+=cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv) # 10  coriolis


                if wind and kchoice==kchoicelist[0]:
                    print("We should not be HERE right now!!!")

                # temporal loop

            print('e_rate =  ',e_rate)
            print('e_hadv =  ',e_hadv)
            print('e_pgr  =  ',e_hdpr)
            print('e_vadv =  ',e_zadv)
            print('e_vdif =  ',e_zdis)
            print('e_cor  =  ',e_cor)
            print('e_wb = ',e_wb)
            print('e_res  =  ',e_res)
            print("end of time loop: ",kchoice)  
        # vertical loop

        print('tmpamp0 ',tmpamp0.mean())
        print('tmpamp1 ',tmpamp1.mean())
        print('tmpamp2 ',tmpamp2.mean())
        print('tmpamp3 ',tmpamp3.mean())
        print('tmpamp4 ',tmpamp4.mean())
        print('tmpamp5 ',tmpamp5.mean())
        print('tmpamp6 ',tmpamp6.mean())
        print('tmpamp7 ',tmpamp7.mean())
        print('tmpamp8 ',tmpamp8.mean())
        print('tmpamp9 ',tmpamp9.mean())
        print('tmpamp10 ',tmpamp10.mean())
   
        tmpamp0=np.fft.fftshift(tmpamp0)
        tmpamp1=np.fft.fftshift(tmpamp1) 
        tmpamp2=np.fft.fftshift(tmpamp2)
        tmpamp3=np.fft.fftshift(tmpamp3) 
        tmpamp4=np.fft.fftshift(tmpamp4)
        tmpamp5=np.fft.fftshift(tmpamp5) 
        tmpamp6=np.fft.fftshift(tmpamp6)
        tmpamp7=np.fft.fftshift(tmpamp7) 
        tmpamp8=np.fft.fftshift(tmpamp8)
        tmpamp9=np.fft.fftshift(tmpamp9) 
        tmpamp10=np.fft.fftshift(tmpamp10)

        #
        ## get 1D spectra
        print('Sum de tmpamp0 (TKE):',tmpamp0.sum())
        [amp0,count,ktmp,dk]=integ_fft2d(tmpamp0,dx)
        [amp1,count,ktmp,dk]=integ_fft2d(tmpamp1,dx)
        [amp2,count,ktmp,dk]=integ_fft2d(tmpamp2,dx)
        [amp3,count,ktmp,dk]=integ_fft2d(tmpamp3,dx)
        [amp4,count,ktmp,dk]=integ_fft2d(tmpamp4,dx)
        [amp5,count,ktmp,dk]=integ_fft2d(tmpamp5,dx)
        [amp6,count,ktmp,dk]=integ_fft2d(tmpamp6,dx)
        [amp7,count,ktmp,dk]=integ_fft2d(tmpamp7,dx)
        [amp8,count,ktmp,dk]=integ_fft2d(tmpamp8,dx)
        [amp9,count,ktmp,dk]=integ_fft2d(tmpamp9,dx)
        [amp10,count,ktmp,dk]=integ_fft2d(tmpamp10,dx)

        if save:
            BDATA=[amp0,amp1,amp2,amp3,amp4,amp5,amp6,amp7,amp8,amp9,amp10,ktmp,count,kchoicelist]
            file1=open(filename,'w')
            pk.dump(BDATA,file1)
            file1.close()
#            ADVECTION=[C4_hadvection_u,C4_hadvection_v,UP3_hadvection_u,UP3_hadvection_u]
#            file2=open("check_diffusion.pck",'w')
#             pk.dump(ADVECTION,file2)
#            file2.close()

    if loaded:
        file1=open(filename,'r')
        BDATALOADED=pk.load(file1)
        file1.close()
        [amp0,amp1,amp2,amp3,amp4,amp5,amp6,amp7,amp8,amp9,amp10,ktmp,count,kchoicelist]=BDATALOADED
    #
    amp11=amp4-(amp1+amp5+amp6+amp7+amp10)
    istr=0
    print(' amp0 rebuilt =  ',amp0.mean())
    #print(' amp12 rebuilt =  ',sum(amp12[istr:]))
    print('e_rate rebuilt =  ',amp4.mean())
    print('e_hdif rebuilt =  ',amp3.mean())
    print('e_hadvup3 rebuilt =  ',amp2.mean())
    print('e_hadc4 rebuilt =  ',amp1.mean())
    print('e_pgr  rebuilt =  ',amp5.mean())
    print('e_vadv rebuilt =  ',amp7.mean())
    print('e_vdif rebuilt =  ',amp6.mean())
    print('e_cor  rebuilt =  ',amp10.mean())
    print('e_res  rebuilt =  ',amp11.mean())
    print('baroclinic conversion rebuilt = ',amp9.mean())
    #
    # #############################################################################
    # #
    # # PLOT
    # #
    #

    dk=np.zeros(ktmp[istr:].shape)
    dk[:-1]=ktmp[istr+1:]-ktmp[istr:-1]
    dk[-1]=dk[-2]
    K=ktmp[istr:]
    A=[]
    
    cff_scale=1.0
    number1=cff_scale

    from math import log10 
    expon=-int(np.floor(log10(abs(number1))))
    coef=number1*10.**expon

    A+=[cff_scale*amp1[istr:]/dk[:]]
    A+=[cff_scale*amp2[istr:]/dk[:]]
    A+=[cff_scale*amp3[istr:]/dk[:]]
    A+=[cff_scale*amp4[istr:]/dk[:]]
    A+=[cff_scale*amp5[istr:]/dk[:]]
    A+=[cff_scale*amp6[istr:]/dk[:]]
    A+=[cff_scale*amp7[istr:]/dk[:]]
    A+=[cff_scale*amp8[istr:]/dk[:]]
    A+=[cff_scale*amp9[istr:]/dk[:]]
    A+=[cff_scale*amp10[istr:]/dk[:]]
    A+=[cff_scale*amp11[istr:]/dk[:]]
    #if wind:
    #  A+=cff_scale*amp22[istr:]/dk
    #
    '''Budget plot'''
    fig=p.figure()
    a=fig.add_subplot(1,1,1)
 #   a.plot(K,A[2],'b',label='horizontal diffusion')
    a.plot(K,A[1],'g',label='horizontal advection')
#    a.plot(K,A[0],'g-*',label='horizontal advection C4')
    a.plot(K,A[3],'r',label=r'$\frac{\partial{u}}{\partial{t}}$')
    a.plot(K,A[4]-A[9],'c',label='P3D')

    a.plot(K,A[5],'b--',label='vertical diffusion',)
    a.plot(K,A[6],'g--',label='vertical advection')
    a.plot(K,A[8],'y--',label='buoyancy flux')
    a.plot(K,A[9],'r--',label='Coriolis')
    #a.plot(K,A[7],'c--',label='Residual')
    #a.plot(K,A[10],'k-*',label='Res(spectral)')
    #a.set_yscale('log')
    a.set_xscale('log')
    a.legend()

    pic1name=resol+'budget'+suffixe
    fig.suptitle(resol+' spectral budget')
    fig.savefig(pic1name)

    ''' Plot residual'''
    fig2=p.figure()
    a2=fig2.add_subplot(1,1,1)
    a2.plot(K,A[7],'c--',label='Residual')
    a2.plot(K,A[10],'k-*',label='Res(spectral)')
    #a.set_yscale('log')
    a2.set_xscale('log')
    a2.legend()
    a2.set_xlabel('k [rad/m]')
    a2.set_ylabel('KE Tendencies [$m^3/s^3$]')
    pic2name=resol+'residual'+suffixe
    fig2.suptitle(resol+' residual')
    fig2.savefig(pic2name)


    '''Plot KE slope'''
    fig3=p.figure()
    a3=fig3.add_subplot(1,1,1)
    a3.set_xlabel('k [rad/m]')
    a3.set_ylabel('KE Tendencies [$m^3/s^3$]')
    #a.set_ylim(ymin,ymax)
    #a3.set_xlim(xmin,xmax)
    pic3name=resol+'slope'+suffixe
    fig3.suptitle(resol)
    cff_scale=1e7
    a3.plot(K,amp0/dk,'go-',label='snapshot')
    a3.plot(K,K**(-2)/cff_scale,'r',label=r'$k^{-2}$')
    a3.plot(K,K**(-5./3)/cff_scale,'b--',label=r'$k^{-5/3}$')
    a3.plot(K,K**(-3)/cff_scale/1e3,'c+',label=r'$k^{-3}$')
    a3.set_yscale('log')
    a3.set_xscale('log')
    a3.set_xlabel('Wavenumber k [rad/m]')
    a3.set_title('KE spectrum')
    a3.legend()
    fig3.savefig(pic3name)
    #p.show()



if __name__ == "__main__":
    import sys
    ''' should use option parser'''
    args=sys.argv[1:]
    resol,depth,tlength=args
    compute_budget(resol,float(depth),int(tlength))
    print("That \'s all folks!")
