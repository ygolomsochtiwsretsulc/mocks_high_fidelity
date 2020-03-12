#!/bin/bash

topcat -stilts plot2plane \
   xpix=500 ypix=400 \
   xlog=true ylog=true xlabel='M_{500c} / M_\odot' \
    ylabel='F_X / \; erg \; cm^{-2} \; s^{-1}' fontsize=16 \
    fontweight=bold texttype=latex \
   xmin=1.0E13 xmax=3.0E15 ymin=5.0E-18 ymax=2.0E-11 \
   auxmap=rainbow auxmin=-2 auxmax=1 \
   auxvisible=true auxlabel='\log_{10}(z)' \
   legend=false \
   layer_1=Mark \
      in_1=$MD10/MD10_eRO_CLU.fit \
      x_1=HALO_M500c y_1=CLUSTER_FX_soft aux_1='log10(redshift_R)' \
      shading_1=aux \
      leglabel_1='4: All' \
   layer_2=Function \
      fexpr_2=6e-14 \
      leglabel_2='6e-14 (-13.2)' \
omode=out out=$GIT_AGN_MOCK/figures/MD10/clusters/FX-M500-z.png

topcat -stilts plot2plane \
   xpix=500 ypix=400 \
   xlog=true xlabel='F_X / \; erg \; cm^{-2} \; s^{-1}' \
    ylabel='L_X / \; erg \; s^{-1}' texttype=latex fontsize=16 fontweight=bold \
   xmin=5.0E-18 xmax=2E-11 ymin=41.5 ymax=45.5 \
   auxmap=rainbow auxmin=0.111 auxmax=1.256 \
   auxvisible=true auxlabel='\log_{10}(kT / \; keV)' \
   legend=false \
   layer=Mark \
      in=$MD10/MD10_eRO_CLU.fit \
      x=CLUSTER_FX_soft y=CLUSTER_LX_soft_RF aux='log10(CLUSTER_kT)' \
      shading=aux size=2 \
omode=out out=$GIT_AGN_MOCK/figures/MD10/clusters/LX-FX-kT.png
      
topcat -stilts plot2plane \
   xpix=500 ypix=400 \
   xlog=true xlabel='F_X / \; erg \; cm^{-2} \; s^{-1}' \
    ylabel='L_X / \; erg \; s^{-1}' texttype=latex fontsize=16 fontweight=bold \
   xmin=5.0E-18 xmax=2E-11 ymin=41.1 ymax=45.5 \
   auxmap=rainbow auxmin=0.0 auxmax=1.5 \
   auxvisible=true auxlabel='z' \
   legend=false \
   layer=Mark \
      in=$MD10/MD10_eRO_CLU.fit \
      x=CLUSTER_FX_soft y=CLUSTER_LX_soft_RF aux=redshift_R \
      shading=aux size=2 \
omode=out out=$GIT_AGN_MOCK/figures/MD10/clusters/LX-FX-z.png

topcat -stilts plot2plane \
   xpix=715 ypix=581 \
   xlog=true xlabel='kT / keV' ylabel=CBP_EM0 \
   xmin=1.22 xmax=19.02 ymin=3.54 ymax=7.51 \
   legend=false \
   in=$MD10/MD10_eRO_CLU.fit x=CLUSTER_kT \
    y=CBP_EM0 \
   layer_1=Mark \
      shading_1=auto size_1=2 \
   layer_2=LinearFit \
      color_2=black thick_2=3 \
omode=out out=$GIT_AGN_MOCK/figures/MD10/clusters/EM0-kT-Xoff.png

topcat -stilts plot2plane \
   xpix=715 ypix=581 \
   xlog=true xlabel=HALO_M500c ylabel=CBP_EM0 \
   xmin=9.3E13 xmax=3.566E15 ymin=3.54 ymax=7.51 \
   legend=false \
   in=$MD10/MD10_eRO_CLU.fit x=HALO_M500c \
    y=CBP_EM0 \
   layer_1=Mark \
      shading_1=auto \
   layer_2=LinearFit \
omode=out out=$GIT_AGN_MOCK/figures/MD10/clusters/EM0-mass.png
   
topcat -stilts plot2plane \
   xpix=600 ypix=500 \
   ylog=true xlabel=z ylabel='|\Delta z|/(1+z)' grid=true texttype=latex \
    fontsize=16 \
   xmin=0.005 xmax=1.5 ymin=1.0E-4 ymax=0.07 \
   auxmap=sron auxmin=376 auxmax=12800 \
   auxvisible=true auxlabel='density' \
   legend=false \
   layer=Grid \
      in=$MD10/MD10_eRO_CLU.fit \
      x=redshift_R y='abs((CBP_redshift-redshift_R)/(1+redshift_R))' \
      xbinsize=-101 ybinsize=-93 combine=count-per-unit \
omode=out out=$GIT_AGN_MOCK/figures/MD10/clusters/deltaz-z.png


topcat -stilts plot2plane \
   xpix=600 ypix=500 \
   xlog=true ylog=true xlabel='M_{500c}' \
    ylabel='|\Delta M_{500c}|/M_{500c} - 0.2' grid=true texttype=latex \
    fontsize=16 \
   xmin=1.0E13 xmax=2.0E15 ymin=1.0E-5 ymax=0.1 \
   auxmap=sron auxmin=0.9 auxmax=159.4 \
   auxvisible=true auxlabel='density' \
   legend=false \
   layer=Grid \
      in=$MD10/MD10_eRO_CLU.fit \
      x=HALO_M500c \
       y='abs((pow(10,CBP_M500c)-HALO_M500c)/HALO_M500c)-0.2' \
      xbinsize=-101 ybinsize=-93 combine=count-per-unit \
omode=out out=$GIT_AGN_MOCK/figures/MD10/clusters/deltaM-M.png

topcat -stilts plot2plane \
   xpix=600 ypix=500 \
   xlog=true ylog=true xlabel='X_{off}' ylabel='EM(0)' grid=true \
    texttype=latex fontsize=16 \
   xmin=0.001 xmax=0.532 ymin=3.0E-8 ymax=2.9E-4 \
   auxmap=sron auxmin=1 auxmax=3009 \
   auxvisible=true auxlabel='density' \
   legend=false \
   layer=Grid \
      in=$MD10/MD10_eRO_CLU.fit \
      x=HALO_Xoff/Halo_rvir y='pow(10,-CBP_EM0)' \
      xbinsize=-101 ybinsize=-93 combine=count-per-unit  \
omode=out out=$GIT_AGN_MOCK/figures/MD10/clusters/xoff-emo.png
