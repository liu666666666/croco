###################################################################
#  this program checks the consistence of diagnostics terms using
# /home3/fastnet/pklein/ES_JUNE2006/s1_2km_100l_fpos_525d_t5_diags

###################################################################
import pyroms
import numpy as np
import pickle as pk
import pylab as p

from utilities import gridinfo,rho2d,cut_var,partialx,partialy,tukeywin,integ_fft2d,integ_fft1d
    
import pycomodo as pc
import pycomodo.operators.simple_ops as op
import pycomodo.util.variables as puv


''' Setting main parameters of for the calculation:
    * if mean =1 , we remove the mean of eaxh variable and use the fluctuation instead u'=u-ubar 
    * window =1 we use windowing (useful for non-periodic solution)
    * cff_tukey=1 if windowing is used uses the Tukey method
    * salinity=1 if we uses salinity in density calculation
    * wind=1 if wind forcing is used'''
mmean=0
window=0
cff_tukey=1
salinity=0
wind=0


''' Path of the output files, resolution, box set, vertical level list, iterations... '''
resol='20K'
dx=20.e3

root='/data/models/JET/Last3Months/'+resol
model='jet_last3months'


kchoicelist=range(30)[8:]#[8]
#kchoicelist=[8]
lims=[1,100,1,400]
#it1=639
it1=1
it2=92


'''physical constants and plotting parameters'''
g=9.81
rho0=1024
cff_scale=1.e7
xmin=0.9e-5
xmax=1.e-3
ymin=-1
ymax=1  # limits for plots
print_and_keep=0
his=-10
dia=-9

#loadmode = input('loadmode? (0/no 1/yes) ')


#
loaded=False
save=True
#
#loaded=True
#save=False

filename=resol+'JET.pck'


''' Different terms are already computed, online=1, or needs to be calculated'''
online_T=1
online_P=1
online_Cor=1
online_Dv=1

kchoice=kchoicelist[0]

#we need even dimensions so:

if np.mod(lims[1]-lims[0],2)==1 :
    lims[0]=lims[0]-1

if np.mod(lims[3]-lims[2],2)==1 :
    lims[2]=lims[2]-1

# Lateral grid
xlim=np.arange(lims[0],lims[1])
ylim=np.arange(lims[2],lims[3])
#overload grid_info function

if not loaded:
    
    grd=gridinfo()
    grd.setId(model)
    
    nc=pc.Archive(grd.grdfile)
    
    #nc=pc.netcdf.netcdf_file(grd.grdfile)
    
    
    pm=nc.variables.get('pm')
    pn=nc.variables.get('pn')
    f=nc.variables.get('f')
    hc=grd.hc
    thetas=grd.thetas
    thetab=grd.thetab
    N=grd.N
    method=grd.Method
    
    h0=nc.variables.get('h')
    
    
    
    #    pm=cut_var(pm.data.transpose(),lims)
    #    pn=cut_var(pn.data.transpose(),lims)
    #    #lon=grd.lonrlon=cut_var(lon,lims)
    #    #lat=grd.latrlat=cut_var(lat,lims)
    #    f=cut_var(f.data.transpose(),lims)
    # Vertical grid
    pm=pm[ylim,xlim].T
    pn=pn[ylim,xlim].T
    f=f[ylim,xlim].T
    
    
    #h=cut_var(h0.data,lims)
    h=h0[0,0]
    
    from Preprocessing_tools import zlevs
    
    zr=zlevs.zlevs(h,0,thetas,thetab,hc,N,'r')
    zw=zlevs.zlevs(h,0,thetas,thetab,hc,N,'w')
    
    kk=kchoicelist
    ll=len(kk)
    dz=np.zeros(zw.shape)
    dz[1:]=zw[1:]-zw[:-1]
    dz[0]=dz[1]
    dz_klist=dz[kk]
    dz=dz[kchoice:N]
    #
    kmin=kchoice-1
    kmax=min(N,kchoice+1)
    [Lmin,Lmax,Mmin,Mmax]=lims
    
    
    
    zr=zlevs.zlevs(h0[:],0,thetas,thetab,hc,N,'r')
    zw=zlevs.zlevs(h0[:],0,thetas,thetab,hc,N,'w')
    
    zr=zr.transpose((1,2,0))
    zw=zw.transpose((1,2,0))
    
    dzw=np.zeros(zw.shape)
    dzw[:,:,1:-1]=zr[:,:,1:]-zr[:,:,:-1]
    dzw[:,:,0]=dzw[:,:,1]
    dzw=dzw[Lmin:Lmax,Mmin:Mmax,kmin:kmax]
    
    dzr=np.zeros(zr.shape)
    dzr[:,:,1:]=zw[:,:,2:]-zw[:,:,1:-1]
    dzr[:,:,0]=dzr[:,:,1]
    dzr=dzr[Lmin:Lmax,Mmin:Mmax,kmin:kmax]
    
    
    historyfilename=root+model+'_his.nc'
    diagfilename =root+model+'_diaM.nc'
    
    print('reading netcdf file: %s',historyfilename)    
    #nchis=pyroms.io.Dataset(historyfilename)
    print('reading netcdf file: %s',diagfilename) 
    ncdiags=pyroms.io.Dataset(diagfilename)
    nchis=pc.Archive(historyfilename)
    #ncdiags=pc.Archive(diagfilename)
    
    
    tmpamp0=np.array(0)
    tmpamp1=np.array(0) 
    tmpamp2=np.array(0)
    tmpamp3=np.array(0) 
    tmpamp4=np.array(0)
    tmpamp5=np.array(0)
    tmpamp6=np.array(0)
    tmpamp7=np.array(0)
    tmpamp8=np.array(0)
    tmpamp9=np.array(0)
    tmpamp10=np.array(0)

   
    ik=-1

    U=nchis.variables['u']
    V=nchis.variables['v']
         
    U=op.interpolate(U,'T')
    V=op.interpolate(V,'T')    
    
    for kchoice in kchoicelist:
        ik=ik+1
        #
        # TEMPORAL LOOP ================================
        #
        for it in range(it1-1,it2):
#            if his==-10:
#                ithis=it+1
#            else:
            ithis=it
            
            print('Treating rec #d',it)
            timehis=nchis.variables['scrum_time'][ithis]
            timediaM=ncdiags.variables['scrum_time'][it]
            print('which, in his file, correspond to %d:', int(timehis))
            print('which, in diaM file, correspond to %d:', int(timediaM))
            
            levelhis=nchis.variables['s_rho'][kchoice]
            leveldiaM=ncdiags.variables['s_rho'][kchoice]
            print('which, in his file, correspond to %lf:', levelhis)
            print('which, in diaM file, correspond to %lf:', leveldiaM)
            
            #### CENTERED DIFFERENCE EVALUATION OF THE NL TERMS ####
            #### computed as if in the rhs
#            u=nchis.variables['u'][ithis,kchoice,:,:]
#            #u=rnt_loadvar_partialz(ctlhis,ithis,'u',kchoice,kchoice) # v is 3d
#            v=nchis.variables['v'][ithis,kchoice,:,:]
            
            ##This is ugly...masked array can't be broadcast? u[iths,kchoice,ylim,xlim]
            u=U[ithis,kchoice,:ylim[-1]+1,:xlim[-1]+1].T
            v=V[ithis,kchoice,:ylim[-1]+1,:xlim[-1]+1].T
#            u2rho=rho2d(u)
#            u=cut_var(u2rho,lims)
#            v2rho=rho2d(v,naxis=0)
#            v=cut_var(v2rho,lims)
#            print('shape 1 (%d,%d)',u.shape)
#            print('u[0,0] %lf',u[0,0])
            #advu1=-u*partialx(u,pm)-v*partialy(u,pn)
            #advv1=-u*partialx(v,pm)-v*partialy(v,pn)
    
    #### UPSTREAM BIASED NL TERMS ####
    #### computed as if in the rhs (- sign already included).
    
    #==============================================================================
            advu2=ncdiags.variables['u_xadv'][it,kchoice,:,:]
            advu2=advu2+ncdiags.variables['u_yadv'][it,kchoice,:,:]
            advu2=rho2d(advu2)
            advu2=cut_var(advu2,lims)
    
            advv2=ncdiags.variables['v_xadv'][it,kchoice,:,:]
            advv2=advv2+ncdiags.variables['v_yadv'][it,kchoice,:,:]
            advv2=rho2d(advv2,naxis=0)
            advv2=cut_var(advv2,lims)
            
            
            advu1=advu2
	    advv1=advv2 
    #
    # #### implicit dissipation
    # #### computed as in the rhs
            #dissu=advu2-advu1
            #dissv=advv2-advv1
    #
    # #### explicit dissipation
    # #### computed as if in the rhs (- sign already included).
            dissu2=ncdiags.variables['u_hmix'][it,kchoice,:,:]
            dissu2=rho2d(dissu2)
            dissu2=cut_var(dissu2,lims)
    
            dissv2=ncdiags.variables['v_hmix'][it,kchoice,:,:]
            dissv2=rho2d(dissv2,naxis=0)
            dissv2=cut_var(dissv2,lims)
    
            dissu=dissu2
            dissv=dissv2
    #
    # #==============================================================================
            if wind:
                print("let's not do that right now!!!")
                exit
    
    #
    # #### temporal derivative
            if online_T==1:
                dudt=ncdiags.variables['u_rate'][it,kchoice,:,:]
                dudt=rho2d(dudt)
                dudt=cut_var(dudt,lims)
                dudt2=dudt
                dvdt=ncdiags.variables['v_rate'][it,kchoice,:,:]
                dvdt=rho2d(dvdt,naxis=0)
                dvdt=cut_var(dvdt,lims)
    
            else:
                print("let's not do that right now!!!")
                exit
    
    
    #
    # #### wb term
            cffb=-g/rho0 # conversion from density to buoyancy
            w=nchis.variables['w'][ithis,kchoice,:,:]
            w=cut_var(w.transpose(),lims)
            #print("show us whatis w.mean(): #lf",w.mean())
            w=w-w.mean()
    
            if salinity:
                T=nchis.variables['temp'][ithis,kchoice,:,:]
                S=nchis.variables['salt'][ithis,kchoice,:,:]
                #rho=rho_potential(T,S)
            else:
                rho=nchis.variables['temp'][ithis,kchoice,:,:]
                
    #
            rho=cut_var(rho.transpose(),lims)
            b=cffb*(rho-rho.mean()) # we remove the mean of b because it is not clear how :)...effectivement c'est pas clair Y.S.
            
    
    #
    # #### pressure term
            if online_P==1:
                dxp=ncdiags.variables['u_Prsgrd'][it,kchoice,:,:]
                dxp=rho2d(dxp)
                dxp=cut_var(dxp,lims)
    
                dyp=ncdiags.variables['v_Prsgrd'][it,kchoice,:,:]
                dyp=rho2d(dyp,naxis=0)
                dyp=cut_var(dyp,lims)
            else:
                print("let's not do that right now!!!")
                exit
    #
    # # Coriolis
    # #
            if online_Cor==1:
                ucor=ncdiags.variables['u_cor'][it,kchoice,:,:]
                ucor=rho2d(ucor)
                ucor=cut_var(ucor,lims)
                vcor=ncdiags.variables['v_cor'][it,kchoice,:,:]
                vcor=rho2d(vcor,naxis=0)
                vcor=cut_var(vcor,lims)
            else:
                print("let's not do that right now!!!")
                exit
                #dndx=np.zeros(pm.shape)
                #dmde=np.zeros(pm.shape)
    #     dndx(2:end-1,:)=0.5./pn(3:end,:)-0.5./pn(1:end-2,:)
    #     dmde(:,2:end-1)=0.5./pm(:,3:end)-0.5./pm(:,1:end-2)
    #     dndx(1,:)=dndx(2,:) dndx(end,:)=dndx(end-1,:)
    #     dmde(:,1)=dmde(:,2) dmde(:,end)=dmde(:,end-1)
    #     cff=f+(v.*dndx-u.*dmde).*pm.*pn
    #     ucor=cff.*v
    #     vcor=-cff.*u
    #     #ucor(2:end-1,:)=0.5*ucor(2:end-1,:)+0.25*(ucor(1:end-2,:)+ucor(3:end,:))
    #     #vcor(:,2:end-1)=0.5*vcor(:,2:end-1)+0.25*(vcor(:,1:end-2)+vcor(:,3:end))
    
    #
    # #### vertical dissipation
            if online_Dv==1:
                ukpp=ncdiags.variables['u_vmix'][it,kchoice,:,:]
                ukpp=rho2d(ukpp)
                ukpp=cut_var(ukpp,lims)
    
                vkpp=ncdiags.variables['v_vmix'][it,kchoice,:,:]
                vkpp=rho2d(vkpp,naxis=0)
                vkpp=cut_var(vkpp,lims)
            else:
                print("let's not do that right now!!!")
                exit
    
    #
    # #### vertical advection
            zadvu=ncdiags.variables['u_vadv'][it,kchoice,:,:]
            zadvu=rho2d(zadvu)
            zadvu=cut_var(zadvu,lims)
    
            zadvv=ncdiags.variables['v_vadv'][it,kchoice,:,:]
            zadvv=rho2d(zadvv,naxis=0)
            zadvv=cut_var(zadvv,lims)
    
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
            resu=rho2d(resu)
            resu=cut_var(resu,lims)
    
            resv=-ncdiags.variables['v_rate'][it,kchoice,:,:]+\
                ncdiags.variables['v_xadv'][it,kchoice,:,:]+\
                ncdiags.variables['v_yadv'][it,kchoice,:,:]+\
                ncdiags.variables['v_cor'][it,kchoice,:,:]+\
                ncdiags.variables['v_vmix'][it,kchoice,:,:]+\
                ncdiags.variables['v_hmix'][it,kchoice,:,:]+\
                ncdiags.variables['v_vadv'][it,kchoice,:,:]+\
                ncdiags.variables['v_Prsgrd'][it,kchoice,:,:]
    
            resv=-resv
            resv=rho2d(resv,naxis=0)
            resv=cut_var(resv,lims)
    
            print('shape 2 (%d,%d)',u.shape)
    # ###############################################################
    #
            if it==it1-1 and kchoice==kchoicelist[0]:
                e_rate=0
                e_hdis=0
                e_hadv=0
                e_hdpr=0
                e_zadv=0
                e_zdis=0
                e_cor=0
                e_res=0
                e_wb=0
    
                (L,M)=u.shape # works for non squared grids
                if window==1:
                    cff_tukey=0.25
                    wdw1=tukeywin(L,cff_tukey)
                    wdw2=tukeywin(M,cff_tukey)
                    print(wdw1.shape)
                    print(wdw2.shape)
                    wdw=np.zeros((wdw1.shape[0],wdw2.shape[0]))
                    #wdw=np.dot(wdw1.transpose(),wdw2)
                    wdw=wdw1.reshape(wdw1.shape[0],1)*wdw2.reshape(1,wdw2.shape[0])
    
            print('shape 3 (%d,%d)',u.shape)
            
            if wind and kchoice==kchoicelist[0]:
                print("let's not do that right now!!!")
                exit
    #     # wind work only at the surface
    #     if mmean==1
    #      us=us-mean(mean(us))vs=vs-mean(mean(vs))
    #      sustr=sustr-mean(mean(sustr)) svstr=svstr-mean(mean(svstr))
    #
    #     if window==1
    #      us=us.*wdwvs=vs.*wdw
    #      sustr=sustr.*wdwsvstr=svstr.*wdw
    #
    #
    #   # apply window and remove spatial mean
    #
    # #  e_rate=e_rate+mean(mean(u.*dudt+v.*dvdt))
    # #  e_hadv=e_hadv+mean(mean(u.*advu2+v.*advv2))
    # #  e_hdpr=e_hdpr+mean(mean(u.*dxp  +v.*dyp))
    # #  e_zadv=e_zadv+mean(mean(u.*zadvu+v.*zadvv))
    # #  e_zdis=e_zdis+mean(mean(u.*ukpp+v.*vkpp))
    # #  e_res=e_rate-(e_hadv+e_hdpr+e_zadv+e_zdis)
    #
            if mmean==1:
                print("We should not be HERE right now!!!") # b and w already done
    #     u=u-mean(mean(u))v=v-mean(mean(v))
    #     advu1=advu1-mean(mean(advu1))
    #     advv1=advv1-mean(mean(advv1))
    #     advu2=advu2-mean(mean(advu2))
    #     advv2=advv2-mean(mean(advv2))
    #     dissu=dissu-mean(mean(dissu))
    #     dissv=dissv-mean(mean(dissv))
    #     dudt=dudt-mean(mean(dudt))
    #     dvdt=dvdt-mean(mean(dvdt))
    #     dxp=dxp-mean(mean(dxp))
    #     dyp=dyp-mean(mean(dyp))
    #     ukpp=ukpp-mean(mean(ukpp))
    #     vkpp=vkpp-mean(mean(vkpp))
    #     zadvu=zadvu-mean(mean(zadvu))
    #     zadvv=zadvv-mean(mean(zadvv))
    #     ucor=ucor-mean(mean(ucor))
    #     vcor=vcor-mean(mean(vcor))
    #     resu=resu-mean(mean(resu))
    #     resv=resv-mean(mean(resv))
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
                print('u[0,0] %lf',u[0,0])
    #
    #
            print('shape 4 (%d,%d)',u.shape)
            e_rate=e_rate+np.mean(u*dudt+v*dvdt)
            e_hdis=e_hdis+np.mean(u*dissu+v*dissv)
            e_hadv=e_hadv+np.mean(u*advu2+v*advv2)
            e_hdpr=e_hdpr+np.mean(u*dxp  +v*dyp)
            e_zadv=e_zadv+np.mean(u*zadvu+v*zadvv)
            e_zdis=e_zdis+np.mean(u*ukpp+v*vkpp)
            e_cor=e_cor+np.mean(u*ucor+v*vcor)
            e_wb=e_wb+np.mean(w*b)
            e_res=e_rate-(e_hadv+e_hdpr+e_zadv+e_zdis+e_cor)
    
    #
    #  SPECTRAL PROJECTION
    #
            if it==it1-1 and kchoice==kchoicelist[0]:
                cff=1./((L*M)**2)
                cff1=cff*dz_klist[ik]/sum(dz_klist)
                tmpamp0=np.zeros((L,M))
                tmpamp1=np.zeros((L,M))
                tmpamp2=np.zeros((L,M))
                tmpamp3=np.zeros((L,M))
                tmpamp4=np.zeros((L,M))
                tmpamp5=np.zeros((L,M))
                tmpamp6=np.zeros((L,M))
                tmpamp7=np.zeros((L,M))
                tmpamp8=np.zeros((L,M))
                tmpamp9=np.zeros((L,M))
                tmpamp10=np.zeros((L,M))
                tmpamp22=np.zeros((L,M))
                tmpamp12=np.zeros(L)
                tmpamp13=np.zeros(L)
    
    #
            print('shape 5 (%d,%d)',tmpamp0.max())
            
            fcoefu=np.fft.fft2(u)
            fcoefv=np.fft.fft2(v)
            tmpamp0=tmpamp0+cff1*np.real(np.conj(fcoefu)*fcoefu+np.conj(fcoefv)*fcoefv)        #  0  kinetic energy
    #
            for i in range(M):
                u1D=u[:,i]
                v1D=v[:,i]
                fcoefu1D=np.fft.fft(u1D)
                fcoefv1D=np.fft.fft(v1D)
                coef=(dz_klist[ik]/(sum(dz_klist)*L**2)).flatten()
                tmpamp12=tmpamp12+coef*np.real(np.conj(fcoefu1D)*fcoefu1D+np.conj(fcoefv1D)*fcoefv1D) 
            print('shape tmpamp12 (%d,%d)',tmpamp12.shape)
            tmpamp12=tmpamp12/M
            tmpamp13=tmpamp0.mean(axis=1)
            print('shape tmpamp13 (%d,%d)',tmpamp13.shape)
            
            fcoeftmpu=np.fft.fft2(advu1)
            fcoeftmpv=np.fft.fft2(advv1)
            tmpamp1=tmpamp1+cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  1  C4 horiz advection
    
            fcoeftmpu=np.fft.fft2(advu2)
            fcoeftmpv=np.fft.fft2(advv2)
            tmpamp2=tmpamp2+cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  2  UP3 horiz advection
            
            fcoeftmpu=np.fft.fft2(dissu)
            fcoeftmpv=np.fft.fft2(dissv)
            tmpamp3=tmpamp3+cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  3  horiz dissipation
            
            fcoeftmpu=np.fft.fft2(dudt)
            fcoeftmpv=np.fft.fft2(dvdt)
            tmpamp4=tmpamp4+cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  4  time rate of change
            
            fcoeftmpu=np.fft.fft2(dxp)
            fcoeftmpv=np.fft.fft2(dyp)
            tmpamp5=tmpamp5+cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  5  pressure gradient
            
            fcoeftmpu=np.fft.fft2(ukpp)
            fcoeftmpv=np.fft.fft2(vkpp)
            tmpamp6=tmpamp6+cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  6  vert dissipation
            
            fcoeftmpu=np.fft.fft2(zadvu)
            fcoeftmpv=np.fft.fft2(zadvv)
            tmpamp7=tmpamp7+cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  7  vert advection
            
            fcoeftmpu=np.fft.fft2(resu)
            fcoeftmpv=np.fft.fft2(resv)
            tmpamp8=tmpamp8+cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv)   #  8  residual
            
            fcoeftmpu=np.fft.fft2(w)
            fcoeftmpv=np.fft.fft2(b)
            tmpamp9=tmpamp9+cff1*np.real(np.conj(fcoeftmpu)*fcoeftmpv)                        #  9  baroclinic conversion wb
            
            fcoeftmpu=np.fft.fft2(ucor)
            fcoeftmpv=np.fft.fft2(vcor)
            tmpamp10=tmpamp10+cff1*np.real(np.conj(fcoefu)*fcoeftmpu+np.conj(fcoefv)*fcoeftmpv) # 10  coriolis
    
    
            if wind and kchoice==kchoicelist[0]:
                print("We should not be HERE right now!!!")
    #     # wind work only at the surface
    #     fcoefus=fft2(us)
    #     fcoefvs=fft2(vs)
    #     fcoeftmpu=fft2(sustr)
    #     fcoeftmpv=fft2(svstr)
    #     tmpamp22=tmpamp22+cff*real(conj(fcoefus).*fcoeftmpu+conj(fcoefvs).*fcoeftmpv)
    #     tmpamp22=1/(it2-it1+1)*tmpamp22
    #     tmpamp22=reorganize_fft(tmpamp22)
    #     method=2
    #     [amp22,count,ktmp,dk]=integ_fft2d(tmpamp22,dx,L,M,method)
    #
    #
    # temporal loop
    
        print('e_rate = %lf ',e_rate)
        print('e_hadv = %lf ',e_hadv)
        print('e_pgr  = %lf ',e_hdpr)
        print('e_vadv = %lf ',e_zadv)
        print('e_vdif = %lf ',e_zdis)
        print('e_cor  = %lf ',e_cor)
        print('e_wb =%lf ',e_wb)
        print('e_res  = %lf ',e_res)
        
        print('shape 6 %lf',tmpamp0.sum())
        print('u[0,0] %lf',u[0,0])
    
    
        tmpamp0=1./(it2-it1+1)*tmpamp0
        tmpamp1=1./(it2-it1+1)*tmpamp1 
        tmpamp2=1./(it2-it1+1)*tmpamp2
        tmpamp3=1./(it2-it1+1)*tmpamp3 
        tmpamp4=1./(it2-it1+1)*tmpamp4
        tmpamp5=1./(it2-it1+1)*tmpamp5
        tmpamp6=1./(it2-it1+1)*tmpamp6
        tmpamp7=1./(it2-it1+1)*tmpamp7 
        tmpamp8=1./(it2-it1+1)*tmpamp8
        tmpamp9=1./(it2-it1+1)*tmpamp9
        tmpamp10=1./(it2-it1+1)*tmpamp10
        tmpamp12=1./(it2-it1+1)*tmpamp12
        tmpamp13=1./(it2-it1+1)*tmpamp13      
        
    # vertical loop
    print('tmpamp0 %lf',tmpamp0.mean())
    print('tmpamp1 %lf',tmpamp1.mean())
    print('tmpamp2 %lf',tmpamp2.mean())
    print('tmpamp3 %lf',tmpamp3.mean())
    print('tmpamp4 %lf',tmpamp4.mean())
    print('tmpamp5 %lf',tmpamp5.mean())
    print('tmpamp6 %lf',tmpamp6.mean())
    print('tmpamp7 %lf',tmpamp7.mean())
    print('tmpamp8 %lf',tmpamp8.mean())
    print('tmpamp9 %lf',tmpamp9.mean())
    print('tmpamp10 %lf',tmpamp10.mean())
    print('shape 6.5 %lf',tmpamp0.sum())
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
    tmpamp12=np.fft.fftshift(tmpamp12)
    tmpamp13=np.fft.fftshift(tmpamp13)
    
    print('tmpamp0 %lf',tmpamp0.mean())
    print('tmpamp1 %lf',tmpamp1.mean())
    print('tmpamp2 %lf',tmpamp2.mean())
    print('tmpamp3 %lf',tmpamp3.mean())
    print('tmpamp4 %lf',tmpamp4.mean())
    print('tmpamp5 %lf',tmpamp5.mean())
    print('tmpamp6 %lf',tmpamp6.mean())
    print('tmpamp7 %lf',tmpamp7.mean())
    print('tmpamp7 %lf',tmpamp7.mean())
    print('tmpamp8 %lf',tmpamp8.mean())
    print('tmpamp9 %lf',tmpamp9.mean())
    print('tmpamp10 %lf',tmpamp10.mean())
    print('tmpamp12 %lf',tmpamp12.mean())
    #
    ## get 1D spectra
    method=2
    print('shape 7 %lf',tmpamp0.sum())
    print(' tmpamp0 before rebuilding = %lf ',sum(sum(tmpamp0)))
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
    
    
    #[amp12,count,ktmp,dk]=integ_fft1d(tmpamp12,dx)
    #[amp13,count,ktmp,dk]=integ_fft1d(tmpamp13,dx)
    #==============================================================================
    # eval(['save SPECTRAL_KE_PK_2d_t' num2str(it1) '_to_t' num2str(it2) '.mat amp0 amp1 amp2 amp3 amp4 amp5 amp6 amp7 amp8 amp9 amp10 ktmp count kchoicelist'])
    #
    if save:
        BDATA=[amp0,amp1,amp2,amp3,amp4,amp5,amp6,amp7,amp8,amp9,amp10,ktmp,count,kchoicelist]
        file1=open(filename,'w')
        pk.dump(BDATA,file1)
        file1.close()

if loaded:
    file1=open(filename,'r')
    BDATALOADED=pk.load(file1)
    file1.close()
    [amp0,amp1,amp2,amp3,amp4,amp5,amp6,amp7,amp8,amp9,amp10,ktmp,count,kchoicelist]=BDATALOADED
#
amp11=amp4-(amp2+amp3+amp5+amp6+amp7+amp10)
istr=0
print(' amp0 rebuilt = %lf ',sum(amp0[istr:]))
#print(' amp12 rebuilt = %lf ',sum(amp12[istr:]))
print('e_rate rebuilt = %lf ',sum(amp4[istr:]))
print('e_hdif rebuilt = %lf ',sum(amp3[istr:]))
print('e_hadv rebuilt = %lf ',sum(amp2[istr:]))
print('e_hadv2 rebuilt = %lf ',sum(amp1[istr:]))
print('e_pgr  rebuilt = %lf ',sum(amp5[istr:]))
print('e_vadv rebuilt = %lf ',sum(amp7[istr:]))
print('e_vdif rebuilt = %lf ',sum(amp6[istr:]))
print('e_cor  rebuilt = %lf ',sum(amp10[istr:]))
print('e_res  rebuilt = %lf ',sum(amp11[istr:]))
print('baroclinic conversion rebuilt = %lf ',sum(amp9[istr:]))
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

print('u[0,0] %lf',u[0,0])
K=ktmp[istr:]
LK=K.shape[0]
A=[]
#number1=K[0]
number1=cff_scale

from math import log10 
expon=-int(np.floor(log10(abs(number1))))
coef=number1*10.**expon
print expon,coef
#cff_scale=10.**4



A +=[cff_scale*amp1[istr:]/dk]
A+=[cff_scale*amp2[istr:]/dk]
A+=[cff_scale*amp3[istr:]/dk]
A+=[cff_scale*amp4[istr:]/dk]
A+=[cff_scale*amp5[istr:]/dk]
A+=[cff_scale*amp6[istr:]/dk]
A+=[cff_scale*amp7[istr:]/dk]
A+=[cff_scale*amp8[istr:]/dk]
A+=[cff_scale*amp9[istr:]/dk]
A+=[cff_scale*amp10[istr:]/dk]
A+=[cff_scale*amp11[istr:]/dk]
#if wind:
#  A+=cff_scale*amp22[istr:]/dk
#
#
# istr=1
# K(1:istr-1)=[]
fig=p.figure()
a=fig.add_subplot(1,1,1)
a.plot(K,A[2],'b',label='horizontal diffusion')
a.plot(K,A[1],'g',label='horizontal advection')
#a.plot(K,A[0],'g-*',label='horizontal advection 2')
a.plot(K,A[3],'r',label=r'$\frac{\partial{u}}{\partial{t}}$')
a.plot(K,A[4],'c',label='P')
#
# hold on
#p.figure()
a.plot(K,A[5],'b--',label='vertical diffusion',)
a.plot(K,A[6],'g--',label='vertical advection')
#a.plot(K,A[8],'y--',label='buoyancy flux')
a.plot(K,A[9],'r--',label='Coriolis')
#a.plot(K,A[7],'c--',label='Residual')
#a.plot(K,A[10],'k-*',label='Res(spectral)')
#a.set_yscale('log')
a.set_xscale('log')


a.legend()
pic1name=resol+'budget'
fig.suptitle(resol+' spectral budget')

fig.savefig(pic1name)


fig2=p.figure()
a2=fig2.add_subplot(1,1,1)
a2.plot(K,A[7],'c--',label='Residual')
a2.plot(K,A[10],'k-*',label='Res(spectral)')
#a.set_yscale('log')
a2.set_xscale('log')


a2.legend()

strexpon=str(expon)
a2.set_xlabel('k [rad/m]')
a2.set_ylabel('KE Tendencies [$10^{'+strexpon+'} m^3/s^3$]')

pic2name=resol+'residual'
fig2.suptitle(resol+' residual')
fig2.savefig(pic2name)
# hl2=plot(K,A(6,istr:),K,A(7,istr:),K,A(10,istr:),K,A(8,istr:),K,A(11,istr:))
# hold on
# set(hl1,'Linewidth',2)
# set(hl2,'Linewidth',2,'Linestyle','--')
# leg('D_H','A_H','T','P','D_V','A_V','COR','Res','Res(spectra)')

#cff_scale
# [coef,expon] =strread(strrep(sprintf('#E',number1),'E','#'),'#f##f')
# number2=floor(coef)*10^(expon)
# set(gca,'Yscale','linear','Xscale','log')
# i#set(gca,'Ylim',[-1 1],'Xlim',[number2 4e-4])
# set(gca,'Ylim',[ymin ymax],'Xlim',[xmin xmax])
fig3=p.figure()
a3=fig3.add_subplot(1,1,1)
a3.set_xlabel('k [rad/m]')
# [coef,expon] =strread(strrep(sprintf('#E',cff_scale),'E','#'),'#f##f')

a3.set_ylabel('KE Tendencies [$10^{'+strexpon+'} m^3/s^3$]')
#a.set_ylim(ymin,ymax)
#a3.set_xlim(xmin,xmax)
#a.set(gca,'Ylim',[ymin ymax],'Xlim',[xmin xmax]);
pic3name=resol+'slope'
fig3.suptitle(resol)

a3.plot(K,amp0/dk,'go-',label='snapshot')
#a3.plot(K,amp12[istr:]/dk,'g*-',label='1D')
#a3.plot(K,amp13[istr:]/dk,'g+-',label='')
a3.plot(K,K**(-2)/cff_scale,'r',label=r'$k^{-2}$')
a3.plot(K,K**(-5./3)/cff_scale,'b--',label=r'$k^{-5/3}$')
a3.plot(K,K**(-3)/cff_scale,'c+',label=r'$k^{-3}$')
a3.set_yscale('log')
a3.set_xscale('log')
a3.set_xlabel('Wavenumber k [rad/m]')
a3.set_title('KE spectrum')
a3.legend()
fig3.savefig(pic3name)
#p.show()
# set(gca,'fontvar.dime',16)
#
# hl=line([3e-6 2e-3],[0 0])
# hold off
# set(hl,'Linewidth',1,'Linestyle','-','Color','g')
# #
#x=[5000:-1000:2000,1000:-500:500,400:-100:100,80:-10:10,5]
#LinkTopAxisData((6.28*1.e-3)./x,x,'Wavelength [km]')
#
# if print_and_keep:
#  outname=['KE_budget_spectra_online']
#  warning off
#  eval(['print -painter -depsc2 ',outname,'.eps'])
#  warning on
#  eval(['! convert -quality 100 ',outname,'.eps ',outname,'.jpg'])
#  eval(['! mv ',outname,'.eps ',dirout_EPS])
#==============================================================================
# eval(['! mv ',outname,'.jpg ',dirout_JPG])
# eval(['! rm -f ',outname,'.eps'])
# eval(['! rm -f ',outname,'.jpg'])
#

