#!/bin/bash


export SASS_CALVERS=/net/ssc4/work1/georg/erosita/caldb/data/erosita/calvers.indx.sletv3


images=efeds_v2_image_all_03-22.fits
expimages=efeds_v2_image_all_03-22_exp.fits
masks=detmask.fits
bkgimages=efeds_comb_ero_bkg_03_22.fits
cheese=efeds_comb_ero_cheese_03_22.fits
srcmaps_phot=efeds_comb_ero_src_ph_03_22.fits
srcmaps_shape=efeds_comb_ero_src_sh_03_22.fits


lboxlist=efeds_comb_ero_lboxlist_03_22.fits
mboxlist=efeds_comb_ero_mboxlist_03_22.fits
mllist_phot=efeds_comb_ero_mllist_ph_03_22.fits
mllist_shape=efeds_comb_ero_mllist_sh_03_22.fits

ecfs="1.28205E12"


emins="300"
emaxs="2200"
eminkev="0.3"
emaxkev="2.1"


doexp=0
domask=0
doboxl=1
dobkg=1
dobox=1
doml=0
doshape=1
dophot=1


if [ $domask -eq 1 ]
then

ermask expimage="$expimages" \
       detmask=$masks \
       threshold1=0.01 \
       threshold2=1.0 \
       regionfile_flag=no

fi


if [ $doboxl -eq 1 ]
then
erbox images="$images" \
      boxlist="$lboxlist" \
      expimages="$expimages" \
      detmasks="$masks" \
      emin="$emins" \
      emax="$emaxs" \
      hrdef=" " \
      ecf=1.0 \
      nruns=3 \
      likemin=6. \
      boxsize=4 \
      compress_flag="N" \
      bkgima_flag="N" \
      expima_flag="Y" \
      detmask_flag="Y"
fi


if [ $dobkg -eq 1 ]
then
    erbackmap image="${images}"  \
              expimage="${expimages}" \
              boxlist="${lboxlist}" \
              detmask="${masks}" \
              cheesemask="${cheese}" \
              bkgimage=${bkgimages} \
              idband=1 \
              scut=0.0001 \
              mlmin=0 \
              maxcut=0.5 \
              fitmethod=smooth \
              nsplinenodes=36 \
              degree=2 \
              smoothflag=yes \
              smoothval=15. \
              snr=40.  \
              excesssigma=10000. \
              nfitrun=1 \
              cheesemaskflag='N'
fi

if [ $dobox -eq 1 ]
then

erbox images="$images" \
      boxlist="$mboxlist" \
      expimages="$expimages" \
      detmasks="$masks" \
      bkgimages="$bkgimages" \
      emin="$emins" \
      emax="$emaxs" \
      hrdef=" " \
      ecf="$ecfs"  \
      nruns=2 \
      likemin=4 \
      boxsize=4 \
      compress_flag="N" \
      bkgima_flag="Y" \
      expima_flag="Y" \
      detmask_flag="Y"


fi


if [ $doshape -eq 1 ]
then

time ermldet mllist=$mllist_shape \
        boxlist=$mboxlist \
        images="$images" \
        expimages="$expimages" \
        detmasks="$masks" \
        bkgimages="$bkgimages" \
        emin="$emins" \
        emax="$emaxs" \
        hrdef=" " \
        ecf="$ecfs" \
        likemin=5. \
        extlikemin=12. \
        extlike_slope="0.0 0.0" \
        compress_flag="N" \
        cutrad=15. \
        multrad=15. \
        extmin=2.5 \
        extmax=30.0 \
        bkgima_flag="Y" \
        expima_flag="Y" \
        detmask_flag="Y" \
        extentmodel="beta" \
        thres_flag="N" \
        thres_col="like" \
        thres_val=30. \
        nmaxfit=3 \
        nmulsou=2 \
        fitext_flag=yes \
        srcima_flag=yes \
        srcimages="$srcmaps_shape" \
        shapelet_flag=yes \
        photon_flag=no \
        sixte_flag=yes &


fi

if [ $dophot -eq 1 ]
then
time ermldet mllist=$mllist_phot \
        boxlist=$mboxlist \
        images="$images" \
        expimages="$expimages" \
        detmasks="$masks" \
        bkgimages="$bkgimages" \
        emin="$emins" \
        emax="$emaxs" \
        hrdef=" " \
        ecf="$ecfs" \
        likemin=5. \
        extlikemin=12. \
        extlike_slope="0.0 0.0" \
        compress_flag="N" \
        cutrad=15. \
        multrad=15. \
        extmin=2.5 \
        extmax=30.0 \
        bkgima_flag="Y" \
        expima_flag="Y" \
        detmask_flag="Y" \
        extentmodel="beta" \
        thres_flag="N" \
        thres_col="like" \
        thres_val=30. \
        nmaxfit=3 \
        nmulsou=2 \
        fitext_flag=yes \
        srcima_flag=yes \
        srcimages="$srcmaps_phot" \
        shapelet_flag=yes \
        photon_flag=yes \
        sixte_flag=yes &



fi
