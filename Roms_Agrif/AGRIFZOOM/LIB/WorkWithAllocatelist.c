/******************************************************************************/
/*                                                                            */
/*     CONV (converter) for Agrif (Adaptive Grid Refinement In Fortran)       */
/*                                                                            */
/* Copyright or � or Copr. Laurent Debreu (Laurent.Debreu@imag.fr)            */
/*                        Cyril Mazauric (Cyril.Mazauric@imag.fr)             */
/* This software is governed by the CeCILL-C license under French law and     */
/* abiding by the rules of distribution of free software.  You can  use,      */
/* modify and/ or redistribute the software under the terms of the CeCILL-C   */
/* license as circulated by CEA, CNRS and INRIA at the following URL          */
/* "http://www.cecill.info".                                                  */
/*                                                                            */
/* As a counterpart to the access to the source code and  rights to copy,     */
/* modify and redistribute granted by the license, users are provided only    */
/* with a limited warranty  and the software's author,  the holder of the     */
/* economic rights,  and the successive licensors  have only  limited         */
/* liability.                                                                 */
/*                                                                            */
/* In this respect, the user's attention is drawn to the risks associated     */
/* with loading,  using,  modifying and/or developing or reproducing the      */
/* software by the user in light of its specific status of free software,     */
/* that may mean  that it is complicated to manipulate,  and  that  also      */
/* therefore means  that it is reserved for developers  and  experienced      */
/* professionals having in-depth computer knowledge. Users are therefore      */
/* encouraged to load and test the software's suitability as regards their    */
/* requirements in conditions enabling the security of their systems and/or   */
/* data to be ensured and,  more generally, to use and operate it in the      */
/* same conditions as regards security.                                       */
/*                                                                            */
/* The fact that you are presently reading this means that you have had       */
/* knowledge of the CeCILL-C license and that you accept its terms.           */
/******************************************************************************/
/* version 1.3                                                                */
/******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "decl.h"


/******************************************************************************/
/*                            AddIdentToTheAllocateList_1                     */
/******************************************************************************/
/* Firstpass 1                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void AddIdentToTheAllocateList_1(char *nom)
{
   listallocate *newvar;
   listallocate *parcours;
   int out;

   if ( firstpass == 1 ) 
   {
      if ( !AllocateList )
      {
         newvar = (listallocate *)malloc(sizeof(listallocate));
         strcpy(newvar->nomvar,nom);
         strcpy(newvar->subroutine,subroutinename);
         strcpy(newvar->module,curmodulename);
         newvar->suiv = NULL;
         AllocateList = newvar;
      }
      else
      {
         parcours = AllocateList;
         out = 0 ; 
         while ( parcours->suiv && out == 0 )
         {
            if (  !strcasecmp(parcours->nomvar,nom) &&
                  !strcasecmp(parcours->subroutine,subroutinename) &&
                  !strcasecmp(parcours->module,curmodulename) ) out = 1;
            else
               parcours=parcours->suiv;               
         }
         if ( out == 0 ) 
         {
            if (  !strcasecmp(parcours->nomvar,nom) &&
                  !strcasecmp(parcours->subroutine,subroutinename) &&
                  !strcasecmp(parcours->module,curmodulename) ) out = 1;
            else
            {
               /* add the record                                              */
              newvar = (listallocate *)malloc(sizeof(listallocate));
              strcpy(newvar->nomvar,nom);
              strcpy(newvar->subroutine,subroutinename);
              strcpy(newvar->module,curmodulename);
              newvar->suiv = NULL;
              parcours->suiv = newvar;
            }
         }
      }
   }
}


/******************************************************************************/
/*                            AddNameToTheAllocateList_1                      */
/******************************************************************************/
/* Firstpass 1                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void AddNameToTheAllocateList_1(char *nom,char *nommodule)
{
   listallocate *newvar;
   listallocate *parcours;
   int out;

   if ( firstpass == 1 ) 
   {
      if ( !AllocateList )
      {
         newvar = (listallocate *)malloc(sizeof(listallocate));
         strcpy(newvar->nomvar,nom);
         strcpy(newvar->subroutine," ");
         strcpy(newvar->module,nommodule);
         newvar->suiv = NULL;
         AllocateList = newvar;
      }
      else
      {
         parcours = AllocateList;
         out = 0 ; 
         while ( parcours->suiv && out == 0 )
         {
            if (  !strcasecmp(parcours->nomvar,nom) &&
                  !strcasecmp(parcours->subroutine," ") &&
                  !strcasecmp(parcours->module,nommodule) ) out = 1;
            else
               parcours=parcours->suiv;               
         }
         if ( out == 0 ) 
         {
            if (  !strcasecmp(parcours->nomvar,nom) &&
                  !strcasecmp(parcours->subroutine," ") &&
                  !strcasecmp(parcours->module,nommodule) ) out = 1;
            else
            {
               /* add the record                                              */
              newvar = (listallocate *)malloc(sizeof(listallocate));
              strcpy(newvar->nomvar,nom);
              strcpy(newvar->subroutine," ");
              strcpy(newvar->module,nommodule);
              newvar->suiv = NULL;
              parcours->suiv = newvar;
            }
         }
      }
   }
}


/******************************************************************************/
/*                            IsAllocateInThisSubroutine_0                    */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int IsAllocateInThisSubroutine_0()
{
   listallocate *parcours;
   int out;

   out = 0 ; 
   if ( firstpass == 0 ) 
   {
      parcours = AllocateList;
      while ( parcours && out == 0 )
      {
         if ( !strcasecmp(parcours->subroutine,subroutinename)  )
         {
            out = 1 ;
         }
         else parcours=parcours->suiv;               
      }
   }
   return out;
}

/******************************************************************************/
/*                            IsVarAllocatable_0                              */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int IsVarAllocatable_0(char *ident)
{
   listallocate *parcours;
   int out;

   out = 0 ; 
   if ( firstpass == 0 ) 
   {
      parcours = AllocateList;
      while ( parcours && out == 0 )
      {
         if ( !strcasecmp(parcours->nomvar,ident)  ) out = 1 ;
         else parcours=parcours->suiv;               
      }
   }
   return out;
}


/******************************************************************************/
/*                          varisallocatable_0                                */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int varisallocatable_0(char *ident)
{
   listvar *newvar;
   listallocate *newvaralloc;
   int out;

   out =0;
   if (firstpass == 0 )
   {
         newvaralloc = AllocateList;
         while ( newvaralloc && out == 0 )
         {
            if ( !strcasecmp(ident,newvaralloc->nomvar) )  out = 1;
            else newvaralloc = newvaralloc->suiv;
         }      
   }
   return out;
}
