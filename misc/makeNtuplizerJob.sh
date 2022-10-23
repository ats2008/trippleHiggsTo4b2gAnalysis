declare -a tagList=( "fileList/MiniAODs/bkg_ggZTo2BHHTo2B2G_UL17.fls" \
          "fileList/MiniAODs/bkg_ttHHTo2B2G_UL17.fls" \
          "fileList/MiniAODs/bkg_ttWToQQHTo2G_UL17.fls" \
          "fileList/MiniAODs/bkg_ZToBBHHTo2B2G_UL17.fls" \
          "fileList/MiniAODs/bkg_WToQQHHTo2B2G_UL17.fls" \
           )

for i in "${tagList[@]}" ;
    do
        echo 
        echo $i
        y=`echo $i | sed 's/.*bkg_//g' | sed 's/.fls.*//g'` ;
        echo $i $y

        ./misc/condorJobMaker.py \
            cms_configs/genNtuplize.py \
            $i \
            results/mc/ntuples/$y/v1p3/ \
            1000 \
            2 \
            -1 \
            $y \
            150
    done
#
#
#./misc/condorJobMaker.py \
#        cms_configs/genNtuplize.py \
#        fileList/c3_5_c4_100_HHHto4b2gamma.fls \
#        results/mc/v1p3/c3_5_c4_100_HHHto4b2gamma/ \
#        10000 \
#        1 \
#        -1 \
#        c3_5_c4_100 \
#        250
