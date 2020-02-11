#!/bin/bash
mkdir hypo.output
cd hypo.output/
# this is for all G+ and G-
awk 'BEGIN {OFS="\t"} {print $3}' ../../step2/hypo.output/hypo.G+.link.zscore.txt > hypo.G+.link.geneID.txt
# let's make zscore for unique geneID #
cp hypo.G+.link.geneID.txt hypo.link.geneID.txt
sort hypo.link.geneID.txt | uniq > hypo.link.geneID.uniq.txt
sort -t $'\t' -k 3,3 ../../step2/hypo.output/hypo.zscore.all.genes.rda.txt > hypo.zscore.all.genes.rda.sorted.txt
awk ' NR==FNR { names[$0]++; next } ($3 in names) { file=$3".hypo.zscore.dataframe.fin.txt"; print $0 >>(file); close (file) } ' hypo.link.geneID.uniq.txt hypo.zscore.all.genes.rda.sorted.txt
Rscript ../../scripts/step3/get.permutation.z.score.for.hypoG+.enhancers.bringfiles.realfast.from.step2.R
allrda=$(wc -l ./hypo.G+.output/*perm.txt)
allrdac=$(echo $allrda | awk '{print $1}')
if [[ $allrdac > 0 ]];
rm *fin.txt;
then echo "hypo.G+ permutaion is done";
fi
