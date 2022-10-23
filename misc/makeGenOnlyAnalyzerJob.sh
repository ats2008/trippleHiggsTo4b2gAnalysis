allTags=( \
"kl3_0_kl4_2" \
"kl3_0p0_kl4_0p0" \
"kl3_0p0_kl4_1p0" \
"kl3_1p0_kl4_0p0" \
"kl3_1p0_kl4_100p0" \
"kl3_1p0_kl4_10p0" \
"kl3_1p0_kl4_1p0" \
"kl3_1p0_kl4_20p0" \
"kl3_1p0_kl4_2p0" \
"kl3_1p0_kl4_4p25" \
"kl3_1p0_kl4_500p0" \
"kl3_1p0_kl4_m100p0" \
"kl3_1p0_kl4_m10p0" \
"kl3_1p0_kl4_m20p0" \
"kl3_1p0_kl4_m2p0" \
"kl3_1p0_kl4_m500p0" \
"kl3_2_kl4_0" \
"kl3_2p0_kl4_1p0" \
"kl3_2p0_kl4_3p0" \
"kl3_2p25_kl4_1p0" \
"kl3_2p25_kl4_40" \
"kl3_2p25_kl4_4p25" \
"kl3_2p25_kl4_m40" \
"kl3_3_kl4_0" \
"kl3_3p0_kl4_2p0" \
"kl3_4_kl4_0" \
"kl3_4p0_kl4_5p0" \
"kl3_6p49_kl4_1p0" \
"kl3_6p49_kl4_40p0" \
"kl3_6p49_kl4_m40" \
"kl3_9p0_kl4_10p0" \
"kl3_m0p5_kl4_0p5" \
"kl3_m1p0_kl4_m1p0" \
"kl3_m1p24_kl4_1p0" \
"kl3_m1p24_kl4_40" \
"kl3_m1p24_kl4_m40" \
)


for TAG in ${allTags[@]}; do
    echo Doing $TAG 

    ./misc/condorJobMakerGeneric.py \
           /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/fastSimStudy/CMSSW_10_2_22/genAnalysis/python/genAnalyzer.py \
           fileList/ntuples/genOnly/$TAG.fls \
           misc/runPython.tpl.sh \
           misc/recoAnalyzer.tpl.cfg \
           results/mc/genOnlyAnalysis/v1p4/$TAG/ \
           10000 \
           1 \
           -1 \
           GenOnlyAly_$TAG \
           250
   done 
