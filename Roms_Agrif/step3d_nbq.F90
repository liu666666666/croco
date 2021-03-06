#include "cppdefs.h"
#ifdef NBQ
!
!======================================================================
!                      NBQ-Mode for NH-modeling
!                            Main Routine
!======================================================================
!
!> @note Main web documentation at:
! http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm
!
! DESCRIPTION: 
!
!> @brief SNBQ driver : Non-hydrostatic algorithm with the 
!                       Non-boussinesq solver.
!
!> @details NBQ time step. See SNBQ web pages :
!  http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/Restricted/NH-NBQ/Html_pages
!    Algorithme_NBQ.htm      --> SNBQ algorithm:               
!    Algebrique_SNBQ.htm     --> SNBQ algebric representation
!    Couplage_Numerique.htm  --> Numerical coupling
!    Couplage_Modes_SNBQ.htm --> Coupling:
!    Couplage_Split_NBQ.htm  --> Coupling Splitting
!
! REVISION HISTORY:
!
!> @authors
!> @date 2015 January
!> @todo
!
!======================================================================
!
      subroutine step3d_nbq (Istr,Iend,Jstr,Jend, WORK)

      use module_nh 
      use module_nbq

      implicit none

# include "param_F90.h"
# include "scalars_F90.h"
# include "ocean3d.h"
# include "grid.h"
# include "nbq.h"
      integer Istr,Iend,Jstr,Jend
      real    WORK(PRIVATE_2D_SCRATCH_ARRAY)
      real :: dum_s
# ifndef MPI
      integer mynode,i
      mynode=0
# else
      integer i
      
# endif

# undef DEBUG
!
!-------------------------------------------------------------------
!       Store internal mode density field
!-------------------------------------------------------------------
!
        call density_nbq(4)
!
!-------------------------------------------------------------------
!  Get internal and external forcing terms for nbq equation:
!  ru+rubar (or lambda_ext+lambda_int)
!-------------------------------------------------------------------
!
      call ru_nbq(1)

!-------------------------------------------------------------------
!       Implicit part: system setup
!-------------------------------------------------------------------
!
# ifdef NBQ_IMP
      if (iif.eq.1.and.ifl_imp_nbq.eq.1) call implicit_nbq (1)
# endif
!
!*******************************************************************
!*******************************************************************
!              NBQ mode iteration (main loop)
!*******************************************************************
!*******************************************************************

      do iteration_nbq=1,iteration_nbq_max

# ifdef DEBUG
        print *,'STEP3D_NBQ: it_nbq = ',iteration_nbq
# endif


!.......1st iteration NBQ: treats NBQ
        if (iteration_nbq>1) then
!        Receive
# ifndef NBQ_IMP
         call parallele_nbq(151)
         call parallele_nbq(152)
         call parallele_nbq(153)
# else
         call parallele_nbq(153)
# endif
!        Momentum equation: switch time indices (move forward)
         call ru_nbq(7)
        endif
!
!-------------------------------------------------------------------
!      Compute divergence term (AMUX):
!          qdm_nbq_a ==> div_nbq_a
!-------------------------------------------------------------------
!
        call density_nbq(6)
!
!-------------------------------------------------------------------
!      Compute Second viscosity (product mat*vect):
!            div_nbq_a ==> rhsd2_nbq
!-------------------------------------------------------------------
!
        call viscous_nbq (1)
!
!-------------------------------------------------------------------
!      Message passing , indice switch (iif=1)
!-------------------------------------------------------------------
!
!  Send
          call parallele_nbq(7)
!
!-------------------------------------------------------------------
!      Acoustic wave emission
!-------------------------------------------------------------------
!
!       call density_nbq(11)
!
!-------------------------------------------------------------------
!      Compute pressure gradient and gravity terms (AMUX)
!                rhp_nbq_a  ==> rhs1_nbq
!-------------------------------------------------------------------
!
        call ru_nbq(6)
!
!-------------------------------------------------------------------
!      Horizontal Momentum equation: leapfrog time stepping
!         If explicit: (x,y,z) is dealt with here
!-------------------------------------------------------------------
!
!  XI-Direction:
!
        do l_nbq = nequ_nh(1)+1,nequ_nh(6)
          dum_s =             soundspeed_nbq**2  * rhs1_nbq(l_nbq)     &
                             - visc2_nbq_a(l_nbq) * rhsd2_nbq(l_nbq)   
          qdm_nbq_a(l_nbq,2) = qdm_nbq_a(l_nbq,0)  + 2.*dtnbq*(        &
                               dum_s                                   &
                             + dqdmdt_nbq_a(l_nbq) )  
          rhssum_nbq_a(l_nbq,2) = rhssum_nbq_a(l_nbq,2)  +  dum_s    
        enddo 
!
!  U-momentum open boundary conditions
!
# ifdef OBC_NBQ
!        call unbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
!
!  Message passing: Send U (51) 
!
        call parallele_nbq(51)
!
!  ETA-Direction:
!
        do l_nbq = neqv_nh(1)+1,neqv_nh(6)  
          dum_s =             soundspeed_nbq**2  * rhs1_nbq(l_nbq)     &
                             - visc2_nbq_a(l_nbq) * rhsd2_nbq(l_nbq)   
          qdm_nbq_a(l_nbq,2) = qdm_nbq_a(l_nbq,0)  + 2.*dtnbq*(        &
                               dum_s                                   &
                             + dqdmdt_nbq_a(l_nbq) )  
          rhssum_nbq_a(l_nbq,2) = rhssum_nbq_a(l_nbq,2)  +  dum_s    
        enddo 
!
!  V-momentum open boundary conditions
!
# ifdef OBC_NBQ
!        call vnbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
!
!  Message passing: Send V (52) 
!
         call parallele_nbq(52)
!
!-------------------------------------------------------------------
!      Vertical Momentum equation: leapfrog time stepping
!         If explicit: (x,y,z) is dealt with here
!         If implicit: (x,y)   only
!-------------------------------------------------------------------
!
# ifndef NBQ_IMP
!
!  Z-Direction: Explicit
!
           do l_nbq = neqw_nh(1)+1,neqw_nh(2)
             dum_s =             soundspeed_nbq**2  * rhs1_nbq(l_nbq)  &
                                - visc2_nbq_a(l_nbq) * rhsd2_nbq(l_nbq)
             qdm_nbq_a(l_nbq,2) = qdm_nbq_a(l_nbq,0)  + 2.*dtnbq*(     &
                                  dum_s                                &
                                + dqdmdt_nbq_a(l_nbq) )  
             rhssum_nbq_a(l_nbq,2) = rhssum_nbq_a(l_nbq,2)  +  dum_s
           enddo 
# else
!
!  Z-Direction: Implicit
!
           call parallele_nbq(151)  ! u only 
           call parallele_nbq(152)  ! v only 
           call implicit_nbq (2)
           call implicit_nbq (3)
# endif
!
!      Vertical momentum open boundary conditions
!
# ifdef OBC_NBQ
!        call wnbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
!
!-------------------------------------------------------------------
!      Message passing 
!-------------------------------------------------------------------
!
!  Send
        call parallele_nbq(53)   ! w only 

!  Receive
        call parallele_nbq(17)
!
!-------------------------------------------------------------------
!      Mass equation: leapfrog time stepping:
!-------------------------------------------------------------------
!
#ifndef NBQ_DRHODT
        do l_nbq = 1 , neqcont_nh
          rhp_nbq_a(l_nbq,2) = rhp_nbq_a(l_nbq,0)                  &
                             - div_nbq_a(l_nbq,1) * 2. * dtnbq 
        enddo
#else
        do l_nbq = 1 , neqcont_nh 
          rhp_nbq_a(l_nbq,2) = rhp_nbq_a(l_nbq,0)                  &
                             - div_nbq_a(l_nbq,1) * 2. * dtnbq     &
                             + rhs1_nbq(neqmom_nh(0)+l_nbq)
        enddo
#endif
!
!-------------------------------------------------------------------
!      Density open boundary conditions
!-------------------------------------------------------------------
!
# ifdef OBC_NBQ
!        call rnbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
!
!-------------------------------------------------------------------
!       Mass equation: switch time indices (move forward)
!-------------------------------------------------------------------
!
        call density_nbq(7)
!
!*******************************************************************
!*******************************************************************

      enddo    ! NBQ loop

!*******************************************************************
!*******************************************************************
!
!
!-------------------------------------------------------------------
!......Message passing 
!-------------------------------------------------------------------
!
! Receive
!        call parallele_nbq(15)
# ifndef NBQ_IMP
       call parallele_nbq(151)       
       call parallele_nbq(152)       
       call parallele_nbq(153)  
# else       
       call parallele_nbq(153)  
# endif

!-------------------------------------------------------------------
!......Move forward: momentum
!-------------------------------------------------------------------
!
       call ru_nbq(7)
!
!-------------------------------------------------------------------
!......Set NBQ/EXT coupling terms
!-------------------------------------------------------------------
!
      call ru_nbq(2)
      call density_nbq(2)


      end subroutine step3d_nbq

#else
      subroutine step3d_nbq_empty
      end subroutine step3d_nbq_empty
#endif
