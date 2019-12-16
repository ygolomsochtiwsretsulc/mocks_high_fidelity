#!/bin/bash

java -jar /home/comparat/software/topcat-full.jar -stilts plot2plane \
   xpix=600 ypix=600 \
   ylog=true xlabel='R_MAG / SDSS r band fiber 1.4arcsec magnitude' \
    ylabel='TEXP_D / min' grid=true fontsize=18 \
   xmin=17 xmax=23 ymin=0.5 ymax=5001 \
   auxmin=0.003 auxmax=1.137 \
   auxvisible=true auxlabel=REDSHIFT \
   legend=true  legpos=0.0,1.0 \
   layer_1=Mark \
      in_1=/home/comparat/data/4most/S5/output_catalogue_S5_LR_4May19.fits \
      x_1=R_MAG y_1=TEXP_D aux_1=REDSHIFT \
      shading_1=aux size_1=5 \
   layer_2=Function \
      fexpr_2='pow(10,0.8*x)*6e-17' color_2=blue thick_2=5 dash_2=3,3 \
      leglabel_2='slope 0.8' \
   layer_3=Function \
      fexpr_3='pow(10,0.4*x)*2e-8' color_3=green thick_3=4 \
      leglabel_3='slope 0.4' \
   layer_4=Function \
      fexpr_4=120 color_4=black thick_4=3 dash_4=12,3,3,3 \
      leglabel_4='2 hours' \
   legseq=_2,_3,_4  \
   omode=out out=Mag_TEXP_D.png



java -jar /home/comparat/software/topcat-full.jar -stilts plot2plane \
   xpix=600 ypix=600 \
   ylog=true xlabel='R_MAG / SDSS r band fiber 1.4arcsec magnitude' \
    ylabel='TEXP_D / min' grid=true  fontsize=18 \
   xmin=17 xmax=23 ymin=0.5 ymax=5001 \
   auxmin=0.003 auxmax=1.137 \
   auxvisible=true auxlabel=REDSHIFT \
   legend=true  legpos=0.0,1.0 \
   layer_1=Mark \
      in_1=/home/comparat/data/4most/S5/output_catalogue_S5_LR_4May19.fits \
      x_1=R_MAG y_1=TEXP_G aux_1=REDSHIFT \
      shading_1=aux size_1=5 \
   layer_2=Function \
      fexpr_2='pow(10,0.8*x)*6e-17' color_2=blue thick_2=5 dash_2=3,3 \
      leglabel_2='slope 0.8' \
   layer_3=Function \
      fexpr_3='pow(10,0.4*x)*2e-8' color_3=green thick_3=4 \
      leglabel_3='slope 0.4' \
   layer_4=Function \
      fexpr_4=120 color_4=black thick_4=3 dash_4=12,3,3,3 \
      leglabel_4='2 hours' \
   legseq=_2,_3,_4  \
   omode=out out=Mag_TEXP_G.png

java -jar /home/comparat/software/topcat-full.jar -stilts plot2plane \
   xpix=600 ypix=600 \
   xlabel='REDSHIFT' ylabel='R AB MAG' grid=true fontsize=18 \
   xmin=0.003 xmax=1.249 ymin=14.95 ymax=24.54 \
   legend=true legpos=0.0,1.0 \
   in=/home/comparat/data/4most/S5/output_catalogue_S5_LR_4May19.fits \
    x=REDSHIFT y=R_MAG \
   layer_1=Mark \
      icmd_1='select ID_SUB_SURVEY==2' \
      color_1=blue size_1=5 \
      leglabel_1='4: cluster_GAL' \
   layer_2=Mark \
      icmd_2='select ID_SUB_SURVEY==3' \
      color_2=green size_2=3 \
      leglabel_2='4: filament_GAL' \
   layer_3=Mark \
      icmd_3='select ID_SUB_SURVEY==1' \
      color_3=red size_3=1 \
      leglabel_3='4: BCG' \
   omode=out out=Mag_Redshift.png 
   
java -jar /home/comparat/software/topcat-full.jar -stilts plot2plane \
   ylog=true xlabel='R AB MAG' ylabel='Counts' \
   xmin=14.81 xmax=24.69 ymin=1 ymax=333690 \
   legend=true legpos=0.0,1.0 \
   in=/home/comparat/data/4most/S5/output_catalogue_S5_LR_4May19.fits \
    x=R_MAG \
   layer_1=Histogram \
      icmd_1='select ID_SUB_SURVEY==1' \
      leglabel_1='1: BCG' \
   layer_2=Histogram \
      icmd_2='select ID_SUB_SURVEY==2' \
      color_2=blue \
      leglabel_2='1: cluster_GAL' \
   layer_3=Histogram \
      icmd_3='select ID_SUB_SURVEY==3' \
      color_3=green \
      leglabel_3='1: filament_GAL' \
   omode=out out=Mag_hist.png

