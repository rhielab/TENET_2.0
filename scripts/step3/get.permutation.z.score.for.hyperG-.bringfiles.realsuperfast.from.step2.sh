#!/bin/bash
mkdir hyper.output
cd hyper.output/
# this is for all G+ and G-
awk 'BEGIN {OFS="\t"} {print $3}' ../../step2/hyper.output/hyper.G-.link.zscore.txt > hyper.G-.link.geneID.txt
# let's make zscore for unique geneID #
cp hyper.G-.link.geneID.txt hyper.link.geneID.txt
sort hyper.link.geneID.txt | uniq > hyper.link.geneID.uniq.txt
sort -t $'\t' -k 3,3 ../../step2/hyper.output/hyper.zscore.all.genes.rda.txt > hyper.zscore.all.genes.rda.sorted.txt
awk ' NR==FNR { names[$0]++; next } ($3 in names) { file=$3".hyper.zscore.dataframe.fin.txt"; print $0 >>(file); close (file) } ' hyper.link.geneID.uniq.txt hyper.zscore.all.genes.rda.sorted.txt
Rscript ../../scripts/step3/get.permutation.z.score.for.hyperG-.enhancers.bringfiles.realfast.from.step2.R
allrda=$(wc -l ./hyper.G-.output/*perm.txt)
allrdac=$(echo $allrda | awk '{print $1}')
if [[ $allrdac > 0 ]];
rm *fin.txt;
then echo "hyper.G- permutaion is done";
fi
