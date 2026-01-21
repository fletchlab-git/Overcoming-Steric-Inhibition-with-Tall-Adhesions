#!/bin/sh



dir='../Data/DataFiles/';
# for i in $(ls ${dir} | grep NoBind);do 
#     tmp="$dir$(echo $i | sed 's/Ll_[0-9]*.[0-9]*_//')";
#     num=$(echo $tmp | cut -d '_' -f4)
#     out=$(echo $tmp | sed -e "s/percent_[0-9]*.[0-9]*_/percent_${num}_/")
#     mv $dir$i ${out}
# done

for i in $(ls ${dir} | grep NoBind | grep -v Ll);do 
  outname=$(echo ${dir}$i | sed -e 's/rhoLBT_0.0_/rhoLBT_0.0_Ll_30.0_/' )
  echo ${dir}$i 
  echo $outname
  echo ""
  
  
  mv ${dir}$i $outname
done
