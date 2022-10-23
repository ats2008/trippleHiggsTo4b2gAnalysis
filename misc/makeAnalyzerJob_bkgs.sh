declare -a tagList=(  \
          "fileList/ntuples/bkg_ggZTo2BHHTo2B2G_UL17.fls" \
          "fileList/ntuples/bkg_ttHHTo2B2G_UL17.fls" \
          "fileList/ntuples/bkg_ttWToQQHTo2G_UL17.fls" \
          "fileList/ntuples/bkg_ZToBBHHTo2B2G_UL17.fls" \
          "fileList/ntuples/bkg_WToQQHHTo2B2G_UL17.fls" \
           )

for fls in "${tagList[@]}" ;
    do
        echo 
        echo $fls
        tag=`echo $fls | sed 's/.*bkg_/bkg_/g' | sed 's/.fls.*//g'` ;
        echo $fls $tag

        ./misc/condorJobMakerGeneric.py \
               /grid_mnt/t3storage3/athachay/trippleHiggs/hhhTo4b2g/fastSimStudy/CMSSW_10_6_29/genAnalysis/python/recoAnalyzer.py \
               $fls \
               misc/runPython.tpl.sh \
               misc/recoAnalyzer.tpl.cfg \
               results/mc/analysis/backgrounds/$tag/v1p5/ \
               5000 \
               1 \
               -1 \
               Aly$tag \
               250
        
   done      
