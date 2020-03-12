#!/usr/bin/env bash
ra_cen=136.2
dec_cen=1.5
Version=4
ATTITUDE="./Simput/eFEDS_att_v2.fits"
GTI_FILE="./Simput/eFEDS_att_v2_gti.fits"

run_sixte(){
NAXIS2=$(fkeyprint ${ATTITUDE} NAXIS2 |awk -F '[=/]' '{if(NF>1){gsub("\047","",$2);print $2,$1}}' |awk '{if($2=="NAXIS2")print $1}')
fdump ${ATTITUDE} STDOUT TIME 1,${NAXIS2} prhead=no showcol=no showunit=no showrow=no align=no page=no pagewidth=256 fldsep=" " wrap=no |awk '{if(NF!=0)print $0}'
TSTART=$(fdump ${ATTITUDE} STDOUT TIME 1 prhead=no showcol=no showunit=no showrow=no align=no page=no pagewidth=256 fldsep=" " wrap=no |awk '{if(NF!=0)print $1+1}')
EXPOSURE=$(fdump ${ATTITUDE} STDOUT TIME 1,${NAXIS2} prhead=no showcol=no showunit=no showrow=no align=no page=no pagewidth=256 fldsep=" " wrap=no |awk '{if(NF!=0)print $0}' |awk '{if(NR==1)tstart=$1}END{print $1-tstart-2}')
printf "START:${TSTART} Exp:${EXPOSURE}\n"
[[ ! -d ${SIMU_DIR} ]] && mkdir ${SIMU_DIR}
erosim Simput=${SIMPUT} Prefix="${SIMU_DIR}/${PreName}_" Attitude=${ATTITUDE} RA=${ra_cen} Dec=${dec_cen} TSTART=${TSTART} Exposure=${EXPOSURE} MJDREF=51543.875 dt=1.0 Seed=100 clobber=yes chatter=3 Background="${AddParticleBkg}" GTIFile=${GTI_FILE}
for Nccd in 1 2 3 4 5 6 7
do
  ero_calevents Projection=AIT Attitude=${ATTITUDE} clobber=yes EvtFile=${SIMU_DIR}/${PreName}_ccd${Nccd}_evt.fits eroEvtFile=${SIMU_DIR}/cal_${PreName}_ccd${Nccd}_evt.fits CCDNR=${Nccd} RA=${ra_cen} Dec=${dec_cen}
  fparkey 15 ${SIMU_DIR}/cal_${PreName}_ccd${Nccd}_evt.fits[1] PAT_SEL add=yes
done
}

AddParticleBkg=yes
PreName="gal"
SIMPUT="./Simput/eFEDS_bkg_v2_simput.fits"
SIMU_DIR="./sim_galacticbkg_v${Version}"
run_sixte
SIMPUT="./Simput/eFEDS_bkg_v2_x2_simput.fits"
SIMU_DIR="./sim_galacticbkg2_v${Version}"
run_sixte

AddParticleBkg=no
SIMU_DIR="./sim_cluster_v${Version}"
SIMPUT="./Simput/cwg_clusters/v2/eFEDS_cls_v2_simput.fits"
PreName="cls"
run_sixte

AddParticleBkg=no
SIMU_DIR="./sim_AGN_v${Version}"
SIMPUT="./Simput/eFEDS_AGN_simput_v4.fits"
PreName="AGN"
run_sixte

ls sim_galacticbkg_v${Version}/cal_gal_ccd[1-7]_evt.fits >in.list
[[ $(wc -l in.list|awk '{print $1}') -ne 7 ]] && echo "ERROR" >&2 && exit 1
evtool eventfiles="@in.list" outfile=evt_eFEDS_onlyBkg_v${Version}.fits pattern=15 clobber=yes
radec2xy file=evt_eFEDS_onlyBkg_v${Version}.fits ra0=${ra_cen} dec0=${dec_cen}
mv evt_eFEDS_onlyBkg_v${Version}.fits /data40s/erosim/eFEDS/Sim${Version}_onlyBkg/

ls sim_AGN_v${Version}/cal_AGN_ccd[1-7]_evt.fits sim_galacticbkg_v${Version}/cal_gal_ccd[1-7]_evt.fits >in.list
[[ $(wc -l in.list|awk '{print $1}') -ne 14 ]] && echo "ERROR" >&2 && exit 1
evtool eventfiles="@in.list" outfile=evt_eFEDS_onlyAGN_v${Version}.fits pattern=15 clobber=yes
radec2xy file=evt_eFEDS_onlyAGN_v${Version}.fits ra0=${ra_cen} dec0=${dec_cen}
mv evt_eFEDS_onlyAGN_v${Version}.fits /data40s/erosim/eFEDS/Sim${Version}_onlyAGN/

ls sim_AGN_v${Version}/cal_AGN_ccd[1-7]_evt.fits sim_galacticbkg_v${Version}/cal_gal_ccd[1-7]_evt.fits sim_cluster_v${Version}/cal_cls_ccd[1-7]_evt.fits >in.list
[[ $(wc -l in.list|awk '{print $1}') -ne 21 ]] && echo "ERROR" >&2 && exit 1
evtool eventfiles="@in.list" outfile=evt_eFEDS_v${Version}.fits pattern=15 clobber=yes
radec2xy file=evt_eFEDS_v${Version}.fits ra0=${ra_cen} dec0=${dec_cen}
mv evt_eFEDS_v${Version}.fits /data40s/erosim/eFEDS/Sim${Version}/

ls sim_AGN_v${Version}/AGN_ccd[1-7]_evt.fits >in.list
[[ $(wc -l in.list|awk '{print $1}') -ne 7 ]] && echo "ERROR" >&2 && exit 1
ftmerge "@in.list" rawevt_eFEDS_AGN_v${Version}.fits
cp rawevt_eFEDS_AGN_v${Version}.fits /data40s/erosim/eFEDS/Sim${Version}/
cp rawevt_eFEDS_AGN_v${Version}.fits /data40s/erosim/eFEDS/Sim${Version}_x2/

ls sim_galacticbkg2_v${Version}/cal_gal_ccd[1-7]_evt.fits >in.list
[[ $(wc -l in.list|awk '{print $1}') -ne 7 ]] && echo "ERROR" >&2 && exit 1
evtool eventfiles="@in.list" outfile=evt_eFEDS_onlyBkg_v${Version}_x2.fits pattern=15 clobber=yes
radec2xy file=evt_eFEDS_onlyBkg_v${Version}_x2.fits ra0=${ra_cen} dec0=${dec_cen}
mv evt_eFEDS_onlyBkg_v${Version}_x2.fits /data40s/erosim/eFEDS/Sim${Version}_onlyBkg_x2/

ls sim_AGN_v${Version}/cal_AGN_ccd[1-7]_evt.fits sim_galacticbkg2_v${Version}/cal_gal_ccd[1-7]_evt.fits sim_cluster_v${Version}/cal_cls_ccd[1-7]_evt.fits >in.list
[[ $(wc -l in.list|awk '{print $1}') -ne 21 ]] && echo "ERROR" >&2 && exit 1
evtool eventfiles="@in.list" outfile=evt_eFEDS_v${Version}_x2.fits pattern=15 clobber=yes
radec2xy file=evt_eFEDS_v${Version}_x2.fits ra0=${ra_cen} dec0=${dec_cen}
mv evt_eFEDS_v${Version}_x2.fits /data40s/erosim/eFEDS/Sim${Version}_x2/

cp rawevt_eFEDS_AGN_v${Version}.fits /data40s/erosim/eFEDS/Sim${Version}_x2/

