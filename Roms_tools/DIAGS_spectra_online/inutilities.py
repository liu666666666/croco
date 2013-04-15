# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 17:31:12 2013

@author: soufflet
"""
import numpy as np

def rho_potential(temp,salt):
#
#=======================================================================
#  Copyright (c) 1996 Rutgers University                             ===
#=======================================================================
#                                                                    ===
#  This routine computes density anomaly via equation of state for   ===
#  seawater.                                                         ===
#                                                                    ===
#  On Input:                                                         ===
#                                                                    ===
#     itracer   Switch indicating which potential temperature and    ===
#               salinity to use (integer):                           ===
#                 itracer=0 => Use prognostic variables.             ===
#                 itracer=1 => Use climatology variables.            ===
#                                                                    ===
#  On Output:                                                        ===
#                                                                    ===
#     dena      Potential Density anomaly (kg/m^3).                  ===
#                                                                    ===
#  Reference:                                                        ===
#                                                                    ===
# << This equation of state formulation has been derived by Jackett  ===
#    and McDougall (1992), unpublished manuscript, CSIRO, Australia. ===
#    It computes in-situ density anomaly as a function of potential  ===
#    temperature (Celsius) relative to the surface, salinity (PSU),  ===
#    and depth (meters).  It assumes  no  pressure  variation along  ===
#    geopotential  surfaces,  that  is,  depth  and  pressure  are   ===
#    interchangeable. >>                                             ===
#                                          John Wilkin, 29 July 92   ===
#                                                                    ===
#=======================================================================
    Q0=+999.842594 ; Q1=+6.793952e-2; Q3=-9.095290e-3;
    Q4=+1.001685e-4; Q5=-1.120083e-6; Q6=+6.536332e-9;
    U0=+0.824493   ; U1=-4.08990e-3 ; U2=+7.64380e-5 ;
    U3=-8.24670e-7 ; U4=+5.38750e-9 ; V0=-5.72466e-3 ;
    V1=+1.02270e-4 ; V2=-1.65460e-6 ; W0=+4.8314e-4;

#-----------------------------------------------------------------------
#  Non-linear equation of state, Jackett and McDougall (1992).
#  MODIFIED XAVIER TO REMOVE THE IN SITU CORRECTION - POTENTIAL
#  DENSITY IS COMPUTED HERE
#-----------------------------------------------------------------------
#
#  Compute secant bulk modulus and store into a utility work array.
#  The units are as follows:
#
#  Compute potential density anomaly (kg/m^3).
#
    dena=Q0+temp*(Q1+temp*(Q3+temp*(Q4+temp*(Q5+temp*Q6))))+\
        salt*(U0+temp*(U1+temp*(U2+temp*(U3+temp*U4))))+\
        salt**(3/2)*(V0+temp*(V1+temp*V2))+W0*salt*salt-1000
    return dena
    
    
    

def reorganize_fft(var):
    (L,M)=var.shape
    MM=M/2
    LL=L/2
    newvar=np.zeros(var.shape)
    newvar[:LL,:MM]=var[LL:,MM:]
    newvar[:LL,MM:]=var[LL:,:MM]
    newvar[LL:,:MM]=var[:LL,MM:]
    newvar[LL:,MM:]=var[:LL,:MM]
    return newvar

    