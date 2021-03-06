#include "cppdefs.h"
#if defined NBQ && defined NBQ_IMP
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
!> @details Loops on the NBQ time step. See SNBQ web pages :
! REVISION HISTORY:
!
!> @authors
!> @date 2015 October
!> @todo
!
!======================================================================
!
      subroutine implicit_nbq ( icall )


!      use module_principal  
!      use module_parallele
       use module_nh 
       use module_nbq
!      use module_exp 
# ifdef TRACETXT
      use module_tracetxt_out
# endif
      implicit none

# include "param_F90.h"
# include "scalars_F90.h"
# include "ocean3d.h"
# include "grid.h"
# include "nbq.h"

      integer :: i,j,k
      integer :: icall,flag_i,ierr_i,ifl_ascii,nd_i


!.....Flag pour les sorties ASCII liees aux calculs implicites:
      ifl_ascii = 0


      if (icall.eq.0) then
!*******************************************************************************
! Initialisations
!*******************************************************************************

      qdmimp_nbq  = 0. 
      pdv_nbq     = 0.
      puv_nbq     = 0.
      plv_nbq     = 0.
      pdvsave_nbq = 0.
      puvsave_nbq = 0.
      plvsave_nbq = 0.
      pdvint_nbq  = 1.
      puvint_nbq  = 0.
      plvint_nbq  = 0.

      endif ! icall == 0

      if (icall.eq.-1) then
!*******************************************************************************
!  Initialisation II
!*******************************************************************************

        call amubdg(                                               &
             neqmimp_nbq                                           &
            ,neqcimp_nbq                                           &
            ,neqmimp_nbq                                           &
            ,mimpj_nbq                                             &
            ,mimpi_nbq                                             &
            ,cimpj_nbq                                             &
            ,cimpi_nbq                                             &
            ,iw1_nbq                                               &
            ,nnz_nbq                                               &
            ,iw2_nbq       )

       endif

      if (icall.eq.1) then
!*******************************************************************************
!   Construction du systeme 
!*******************************************************************************
       flag_i = 1

       ptri_nbq = neqmimp_nbq
       call amub2_tri(                                             &
             neqmimp_nbq                                           &
            ,neqcimp_nbq                                           &
            ,flag_i                                                &
            ,mimpv_nbq                                             &  ! Taille!
            ,mimpj_nbq                                             &
            ,mimpi_nbq                                             &
            ,cimpv_nbq(1:)                                         &
            ,cimpj_nbq                                             &
            ,cimpi_nbq                                             &
            ,pimpj_nbq                                             &
            ,pimpi_nbq                                             &
            ,nnz_nbq                                               &
            ,iw2_nbq                                               &
            ,ierr_i      )

      endif  ! icall == 1

      if (icall.eq.2) then
!*******************************************************************************
!   RHS 
!*******************************************************************************

!.......QDM pour le RHS: terme CONTxQDM
        qdmimp_nbq(1:neqmom_nh(1)+neqmom_nh(2)) =  &
            qdm_nbq_a(1:neqmom_nh(1)+neqmom_nh(2),2) 

!.......Produit MatxVect: CONTxQDM
        call amux(                                                &
              neqcont_nh                                          &
             ,qdmimp_nbq(1)                                       &
             ,rhsimp2_nbq(1)                                      &
             ,contv_nh(1:)                                        &
             ,contj_nh(1:)                                        &
             ,conti_nh(1:)       )                           

       rhsimp2_nbq(1:neqcimp_nbq)=  (rhp_nbq_a(1:neqcimp_nbq,1)      &
                                     -rhp_bq_a(1:neqcimp_nbq))       &
                         - dtnbq * rhsimp2_nbq(1:neqcimp_nbq) 

!.......RHS = MOMxRHS(cont) 
        call amux(                                                   &
               neqmimp_nbq                                           &
              ,rhsimp2_nbq(1:neqcimp_nbq)                            & 
              ,rhsimp_nbq                                            &
              ,mimpv_nbq                                             &
              ,mimpj_nbq                                             &
              ,mimpi_nbq       )  

!.....Sauvegarde du second membre &
!     Ajouts des deux termes de l'equation QDM  du RHS: MOM * RHS(cont) + RHS(mom)
!     Equation de continuite:

      rhssum_nbq_a(neqmom_nh(1)+neqmom_nh(2)+1:neqmom_nh(0),2) =                 & 
      rhssum_nbq_a(neqmom_nh(1)+neqmom_nh(2)+1:neqmom_nh(0),2)                   &
                                  + rhsimp_nbq(1:neqmimp_nbq)/dtnbq 

      rhsimp_nbq(1:neqmimp_nbq) =                                                &
        qdm_nbq_a(neqmom_nh(1)+neqmom_nh(2)+1:neqmom_nh(0),1)                    &
                            + rhsimp_nbq(1:neqmimp_nbq)                          &
                 !          - visc2_nbq_a(l_nbq) * rhsd2_nbq(l_nbq) * dtnbq      &  ! Euler instable!
      + dqdmdt_nbq_a(neqmom_nh(1)+neqmom_nh(2)+1:neqmom_nh(0)) * dtnbq
     
      endif     ! icall == 2

      if (icall.eq.3) then
!*******************************************************************************
! Inversion du systeme implicite
!
!      matrice : ( pimpi_nbq, pimpj_nbq, pimpv_nbq)
!
!      rhs     : rhsimp_nbq
!
!      # de lignes: neqmimp_nbq
!
!      # d'elements non nuls: pimpi_nbq(neqimp_nbq+1)-1
!
!
!*******************************************************************************

      if (iif.eq.1.and.iteration_nbq.eq.1) then

!......Sauvegarde de la matrice tri-diagonale
!        necessaire car dgtsv modifie la matrice tri-diag...

       pdvsave_nbq(1:ptri_nbq) = pdv_nbq(1:ptri_nbq)
       plvsave_nbq(1:ptri_nbq) = plv_nbq(1:ptri_nbq)
       puvsave_nbq(1:ptri_nbq) = puv_nbq(1:ptri_nbq)

      else

!......Recopie de la matrice tri-diagonale...

       pdv_nbq(1:ptri_nbq) = pdvsave_nbq(1:ptri_nbq)
       plv_nbq(1:ptri_nbq) = plvsave_nbq(1:ptri_nbq)
       puv_nbq(1:ptri_nbq) = puvsave_nbq(1:ptri_nbq)

      endif
 
!..... Resolution du systeme tri-diagonal:

      nd_i = 1  ! Une seule inversion

      call dgtsv  ( ptri_nbq                  &
                  ,nd_i                       &
                  ,plv_nbq    (1:ptri_nbq-1)  &
                  ,pdv_nbq    (1:ptri_nbq  )  &
                  ,puv_nbq    (1:ptri_nbq-1)  &
                  ,rhsimp_nbq (1:ptri_nbq)    &
                  ,ptri_nbq                   &
                  ,ierr_i )  

!......Recopie de la solution implicite pour w (update) et calcul rhs (sum)

      qdm_nbq_a(neqw_nh(1)+1:neqw_nh(2),2) =  &
                rhsimp_nbq(neqw_nh(1)-neqmom_nh(1)-neqmom_nh(2)+1:neqw_nh(2)-neqmom_nh(1)-neqmom_nh(2))


      endif

 
       end subroutine implicit_nbq 
#else
      subroutine implicit_nbq_empty
      end subroutine implicit_nbq_empty
#endif
