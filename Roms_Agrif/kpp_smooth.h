! $Id$
!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
!
! ROMS_AGRIF website : http://roms.mpl.ird.fr
!======================================================================
!
! PART OF KPP2005 (Shchepetkin et al. 2005)
!
#   ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=J_EXT_RANGE
          hwrk(Istr-1,j)=hwrk(Istr,j)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=J_EXT_RANGE
          hwrk(Iend+1,j)=hwrk(Iend,j)
        enddo
      endif
#   endif
#   ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=I_EXT_RANGE
          hwrk(i,Jstr-1)=hwrk(i,Jstr)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=I_EXT_RANGE
          hwrk(i,Jend+1)=hwrk(i,Jend)
        enddo
      endif
#    ifndef EW_PERIODIC
      if (WESTERN_EDGE.and.SOUTHERN_EDGE) then
        hwrk(Istr-1,Jstr-1)=hwrk(Istr,Jstr)
      endif
      if (WESTERN_EDGE.and.NORTHERN_EDGE) then
        hwrk(Istr-1,Jend+1)=hwrk(Istr,Jend)
      endif
      if (EASTERN_EDGE.and.SOUTHERN_EDGE) then
        hwrk(Iend+1,Jstr-1)=hwrk(Iend,Jstr)
      endif
      if (EASTERN_EDGE.and.NORTHERN_EDGE) then
        hwrk(Iend+1,Jend+1)=hwrk(Iend,Jend)
      endif
#    endif
#   endif

      do j=Jstr,Jend+1
        do i=Istr,Iend+1
          wrk3(i,j)=0.25*(hwrk(i,j)  +hwrk(i-1,j)
     &                   +hwrk(i,j-1)+hwrk(i-1,j-1))
#   ifdef MASKING
     &                   *pmask(i,j)
#   endif
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
#   ifdef MASKING
          cff=0.25*(pmask(i,j)   +pmask(i+1,j)
     &             +pmask(i,j+1) +pmask(i+1,j+1))
#   else
          cff=1.
#   endif
          hwrk(i,j)=(1.-cff)*hwrk(i,j)+
     &              0.25*(wrk3(i,j)  +wrk3(i+1,j)
     &                   +wrk3(i,j+1)+wrk3(i+1,j+1))
#   ifdef MASKING
          hwrk(i,j)=hwrk(i,j)*rmask(i,j)
#   endif
        enddo
      enddo