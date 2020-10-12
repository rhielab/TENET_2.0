#!/bin/bash

cat cg*rda.txt > hypo.zscore.all.genes.rda.txt

Zscore=$(Rscript ../../scripts/step2/get.zscore.from.pval.R)
IFS="\t" <<< "$Zscore"
Zscorecutoff=$(echo $Zscore | awk '{print $2}')

p="-1"
realcutoff=$(echo "$p*$Zscorecutoff" | bc -l)
awk -v var="$realcutoff" 'BEGIN { OFS="\t"} {if ($4 < var) print $1, $2, $3, $4 }' hypo.zscore.all.genes.rda.txt > hypo.G+.link.zscore.txt
realcutoff=$(echo "$Zscorecutoff" | bc -l)
awk -v var="$realcutoff" 'BEGIN { OFS="\t"} {if ($4 > var) print $1, $2, $3, $4 }' hypo.zscore.all.genes.rda.txt > hypo.G-.link.zscore.test.txt
grep -v NaN hypo.G-.link.zscore.test.txt > hypo.G-.link.zscore.cut.txt
grep -v Inf hypo.G-.link.zscore.cut.txt > hypo.G-.link.zscore.txt
rm hypo.G-.link.zscore.cut.txt
rm hypo.G-.link.zscore.test.txt
allrda=$(wc -l hypo.zscore.all.genes.rda.txt)
allrdac=$(echo $allrda | awk '{print $1}')
if [[ $allrdac > 0 ]];
rm cg*rda.txt;
rm ch*rda.txt;
rm rs*rda.txt;
then echo "hypo zscore selection is done";
fi
