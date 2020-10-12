#!/bin/bash
end=$(ls *bed | wc -l)
for ((i=1; i<=end; i++));
do
	thefile=$(ls *bed | head -n $i| tail -n 1)
	intersectBed -a ../../scripts/data/hm450cg_GEO.GR.bed -b $thefile -wa -wb > $thefile.hm450cg_GEO_GR.txt
	awk 'BEGIN { OFS="\t"} {print $4}' $thefile.hm450cg_GEO_GR.txt | sort | uniq > $thefile.probelist.txt
	mv $thefile.probelist.txt ../../scripts/step1/hm450overlap/feature/
	rm *txt
done
cat ../../scripts/step1/hm450overlap/feature/*txt | sort | uniq > ../../scripts/step1/hm450overlap/feature/external.data.feature.hm450.probelist.txt
