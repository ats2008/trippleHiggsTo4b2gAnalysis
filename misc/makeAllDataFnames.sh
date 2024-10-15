 for  i in  results/STORE/eventIDRemade/evtIDXUpdaterData_data_201*  ; 
    do echo $i ;  
        tag=`sed 's#/#_#g'<<<$i`;
        #echo "realpath $i/* > $tag.fls ";
        realpath $i/* > $tag.fls 
        echo "      Made file list : " $tag.fls;
    done
