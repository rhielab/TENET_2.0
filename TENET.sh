#!/bin/bash
echo "Welcome to TENET"
cp settings.txt settings.R
Rscript settings.R
cores=$(head -n 3 settings.txt | tail -n 1)
eval $cores
# all step parameters #
step0=$(head -n 6 settings.txt | tail -n 1)
eval $step0
step1=$(head -n 7 settings.txt | tail -n 1)
eval $step1
step2=$(head -n 8 settings.txt | tail -n 1)
eval $step2
step3=$(head -n 9 settings.txt | tail -n 1)
eval $step3
step4=$(head -n 10 settings.txt | tail -n 1)
eval $step4
step5=$(head -n 11 settings.txt | tail -n 1)
eval $step5
# step 0 #
if [[ $step0 == "F" ]];
then
echo "Thanks for loading files"
fi
if [[ $step0 == "T" ]];
then
echo "step0 started"
cd external.data
cd data
Rscript ../../scripts/step0/make.rdafiles.from.textfiles.R
echo "step0 done"
cd ../../
fi
othercells=$(head -n 24 settings.txt | tail -n 1)
eval $othercells
if [[ $step0 == "T" && $othercells == "T" ]];
then
cd external.data
cd othercells
Rscript ../../scripts/step1/make.rdafiles.from.textfiles.for.othercells.R
cd ../../
fi
# step 1 #
if [[ $step1 == "T" ]];
then
mkdir step1
mkdir './scripts/step1/hm450overlap'
mkdir './scripts/step1/hm450overlap/enhancer'
mkdir './scripts/step1/hm450overlap/NDR'
mkdir './scripts/step1/hm450overlap/feature'
cd ./step1
echo "step1 started"
extENH=$(head -n 36 ../settings.txt | tail -n 1)
eval $extENH
if [[ $extENH == "T" ]];
then 
cd ../external.data/enhancer
bash ../../scripts/step1/make.enhancer.hm450.probelist.then.combine.sh;
cd ../../step1
fi
extNDR=$(head -n 37 ../settings.txt | tail -n 1)
eval $extNDR
if [[ $extNDR == "T" ]];
then
cd ../external.data/NDR
bash ../../scripts/step1/make.NDR.hm450.probelist.then.combine.sh;
cd ../../step1
fi
onlyextFeature=$(head -n 38 ../settings.txt | tail -n 1)
eval $onlyextFeature
if [[ $onlyextFeature == "T" ]];
then 
cd ../external.data/feature
bash ../../scripts/step1/make.feature.hm450.probelist.then.combine.sh;
cd ../../step1
fi
Rscript ../scripts/step1/get.diffmeth.regions.R;
cd ..
echo "step1 done"
fi
# step 2 #
if [[ $step2 == "T" ]];
then
mkdir step2
cd ./step2
echo "step2 started"
findhypo=$(head -n 46 ../settings.txt | tail -n 1)
eval $findhypo
if [[ $findhypo == "T" ]];
then
Rscript ../scripts/step2/get.pairs.z.score.for.hypometh.enhancers.realfast.R
cd ./hypo.output
bash ../../scripts/step2/combine.z.score.hypometh.files.then.cut.sh
cd ..
fi
findhyper=$(head -n 53 ../settings.txt | tail -n 1)
eval $findhyper
if [[ $findhyper == "T" ]];
then
Rscript ../scripts/step2/get.pairs.z.score.for.hypermeth.enhancers.realfast.R
cd ./hyper.output
bash ../../scripts/step2/combine.z.score.hypermeth.files.then.cut.sh
cd ..
fi
cd ..
echo "step2 done"
fi
# step 3 #
if [[ $step3 == "T" ]];
then 
mkdir step3
cd ./step3
echo "step3 started"
findhypoGpos=$(head -n 49 ../settings.txt | tail -n 1)
eval $findhypoGpos
findhypoGneg=$(head -n 50 ../settings.txt | tail -n 1)
eval $findhypoGneg
if [[ $findhypoGpos == "T" && $findhypoGneg == "T" ]];
then
bash ../scripts/step3/get.permutation.z.score.for.hypometh.bringfiles.realsuperfast.from.step2.sh
fi
if [[ $findhypoGpos == "T" && $findhypoGneg == "F" ]];
then
bash ../scripts/step3/get.permutation.z.score.for.hypoG+.bringfiles.realsuperfast.from.step2.sh
fi
if [[ $findhypoGpos == "F" && $findhypoGneg == "T" ]];
then
bash ../scripts/step3/get.permutation.z.score.for.hypoG-.bringfiles.realsuperfast.from.step2.sh
fi
findhyperGpos=$(head -n 56 ../settings.txt | tail -n 1)
eval $findhyperGpos
findhyperGneg=$(head -n 57 ../settings.txt | tail -n 1)
eval $findhyperGneg
if [[ $findhyperGpos == "T" && $findhyperGneg == "T" ]];
then
bash ../scripts/step3/get.permutation.z.score.for.hypermeth.bringfiles.realsuperfast.from.step2.sh
fi
if [[ $findhyperGpos == "T" && $findhyperGneg == "F" ]];
then
bash ../scripts/step3/get.permutation.z.score.for.hyperG+.bringfiles.realsuperfast.from.step2.sh
fi
if [[ $findhyperGpos == "F" && $findhyperGneg == "T" ]];
then
bash ../scripts/step3/get.permutation.z.score.for.hyperG-.bringfiles.realsuperfast.from.step2.sh
fi
cd ..
echo "step3 done"
fi
# step 4 #
if [[ $step3 == "T" && $step4 == "T" ]];
then
mkdir step4
cd ./step4
echo "step4 started"
findhypoGpos=$(head -n 49 ../settings.txt | tail -n 1)
eval $findhypoGpos
if [[ $findhypoGpos == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypometh.G+.pairs.realfast.from.step3.R
fi 
findhypoGneg=$(head -n 50 ../settings.txt | tail -n 1)
eval $findhypoGneg
if [[ $findhypoGneg == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypometh.G-.pairs.realfast.from.step3.R
fi
findhyperGpos=$(head -n 56 ../settings.txt | tail -n 1)
eval $findhyperGpos
if [[ $findhyperGpos == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypermeth.G+.pairs.realfast.from.step3.R
fi
findhyperGneg=$(head -n 57 ../settings.txt | tail -n 1)
eval $findhyperGneg
if [[ $findhyperGneg == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypermeth.G-.pairs.realfast.from.step3.R
fi
cd ..
echo "step4 done"
fi
if [[ $step2 == "T" && $step3 == "F" && $step4 == "T" ]];
then
mkdir step4
cd ./step4
echo "step4 started"
findhypoGpos=$(head -n 49 ../settings.txt | tail -n 1)
eval $findhypoGpos
if [[ $findhypoGpos == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypometh.G+.pairs.realfast.from.step2.R
fi 
findhypoGneg=$(head -n 50 ../settings.txt | tail -n 1)
eval $findhypoGneg
if [[ $findhypoGneg == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypometh.G-.pairs.realfast.from.step2.R
fi
findhyperGpos=$(head -n 56 ../settings.txt | tail -n 1)
eval $findhyperGpos
if [[ $findhyperGpos == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypermeth.G+.pairs.realfast.from.step2.R
fi
findhyperGneg=$(head -n 57 ../settings.txt | tail -n 1)
eval $findhyperGneg
if [[ $findhyperGneg == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypermeth.G-.pairs.realfast.from.step2.R
fi
cd ..
echo "step4 done"
fi
if [[ $step1 == "T" && $step2 == "F" && $step3 == "F" && $step4 == "T" ]];
then
mkdir step4
cd ./step4
echo "step4 started"
findhypoGpos=$(head -n 49 ../settings.txt | tail -n 1)
eval $findhypoGpos
if [[ $findhypoGpos == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypometh.G+.pairs.realfast.from.step1.R
fi 
findhypoGneg=$(head -n 50 ../settings.txt | tail -n 1)
eval $findhypoGneg
if [[ $findhypoGneg == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypometh.G-.pairs.realfast.from.step1.R
fi
findhyperGpos=$(head -n 56 ../settings.txt | tail -n 1)
eval $findhyperGpos
if [[ $findhyperGpos == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypermeth.G+.pairs.realfast.from.step1.R
fi
findhyperGneg=$(head -n 57 ../settings.txt | tail -n 1)
eval $findhyperGneg
if [[ $findhyperGneg == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypermeth.G-.pairs.realfast.from.step1.R
fi
cd ..
echo "step4 done"
fi
if [[ $step1 == "F" && $step2 == "F" && $step3 == "F" && $step4 == "T" ]];
then
mkdir step4
cd ./step4
echo "step4 started"
findhypoGpos=$(head -n 49 ../settings.txt | tail -n 1)
eval $findhypoGpos
if [[ $findhypoGpos == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypometh.G+.pairs.realfast.from.step3.R
fi
findhypoGneg=$(head -n 50 ../settings.txt | tail -n 1)
eval $findhypoGneg
if [[ $findhypoGneg == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypometh.G-.pairs.realfast.from.step3.R
fi
findhyperGpos=$(head -n 56 ../settings.txt | tail -n 1)
eval $findhyperGpos
if [[ $findhyperGpos == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypermeth.G+.pairs.realfast.from.step3.R
fi
findhyperGneg=$(head -n 57 ../settings.txt | tail -n 1)
eval $findhyperGneg
if [[ $findhyperGneg == "T" ]];
then
Rscript ../scripts/step4/get.optimized.hypermeth.G-.pairs.realfast.from.step3.R
fi
cd ..
echo "step4 done"
fi
if [[ $step2 == "T" && $step3 == "T" && $step4 == "F" ]];
then
mkdir step4
cd ./step4
echo "step4 started"
findhypoGpos=$(head -n 49 ../settings.txt | tail -n 1)
eval $findhypoGpos
if [[ $findhypoGpos == "T" ]];
then
mkdir hypo.G+.output
cd hypo.G+.output
ln -s ../../step3/hypo.output/hypo.G+.output/hypo.G+.link.zscore.perm.txt hypo.G+.link.zscore.perm.all.optimized.links.txt 
cd ..
fi
findhypoGneg=$(head -n 50 ../settings.txt | tail -n 1)
eval $findhypoGneg
if [[ $findhypoGneg == "T" ]];
then
mkdir hypo.G-.output
cd hypo.G-.output
ln -s ../../step3/hypo.output/hypo.G-.output/hypo.G-.link.zscore.perm.txt hypo.G-.link.zscore.perm.all.optimized.links.txt
cd ..
fi
findhyperGpos=$(head -n 56 ../settings.txt | tail -n 1)
eval $findhyperGpos
if [[ $findhyperGpos == "T" ]];
then
mkdir hyper.G+.output
cd hyper.G+.output
ln -s ../../step3/hyper.output/hyper.G+.output/hyper.G+.link.zscore.perm.txt hyper.G+.link.zscore.perm.all.optimized.links.txt
cd ..
fi
findhyperGneg=$(head -n 57 ../settings.txt | tail -n 1)
eval $findhyperGneg
if [[ $findhyperGneg == "T" ]];
then
mkdir hyper.G-.output
cd hyper.G-.output
ln -s ../../step3/hyper.output/hyper.G-.output/hyper.G-.link.zscore.perm.txt hyper.G-.link.zscore.perm.all.optimized.links.txt
cd ..
fi
cd ..
echo "step4 done"
fi
if [[ $step2 == "T" && $step3 == "F" && $step4 == "F" ]];
then
mkdir step4
cd ./step4
echo "step4 started"
findhypoGpos=$(head -n 49 ../settings.txt | tail -n 1)
eval $findhypoGpos
if [[ $findhypoGpos == "T" ]];
then
mkdir hypo.G+.output
cd hypo.G+.output
awk 'BEGIN { OFS="\t"} {print $1, $2, $3, $4, "NA"}' ../../step2/hypo.output/hypo.G+.link.zscore.txt > hypo.G+.link.zscore.all.optimized.links.txt
sed -i $'1 i\\\nprobe\tgeneSymbol\tgeneID\tZ.real\tpval' hypo.G+.link.zscore.all.optimized.links.txt
cd ..
fi
findhypoGneg=$(head -n 50 ../settings.txt | tail -n 1)
eval $findhypoGneg
if [[ $findhypoGneg == "T" ]];
then
mkdir hypo.G-.output
cd hypo.G-.output
awk 'BEGIN { OFS="\t"} {print $1, $2, $3, $4, "NA"}' ../../step2/hypo.output/hypo.G-.link.zscore.txt > hypo.G-.link.zscore.all.optimized.links.txt
sed -i $'1 i\\\nprobe\tgeneSymbol\tgeneID\tZ.real\tpval' hypo.G-.link.zscore.all.optimized.links.txt
cd ..
fi
findhyperGpos=$(head -n 56 ../settings.txt | tail -n 1)
eval $findhyperGpos
if [[ $findhyperGpos == "T" ]];
then
mkdir hyper.G+.output
cd hyper.G+.output
awk 'BEGIN { OFS="\t"} {print $1, $2, $3, $4, "NA"}' ../../step2/hyper.output/hyper.G+.link.zscore.txt > hyper.G+.link.zscore.all.optimized.links.txt
sed -i $'1 i\\\nprobe\tgeneSymbol\tgeneID\tZ.real\tpval' hyper.G+.link.zscore.all.optimized.links.txt
cd ..
fi
findhyperGneg=$(head -n 57 ../settings.txt | tail -n 1)
eval $findhyperGneg
if [[ $findhyperGneg == "T" ]];
then
mkdir hyper.G-.output
cd hyper.G-.output
awk 'BEGIN { OFS="\t"} {print $1, $2, $3, $4, "NA"}' ../../step2/hyper.output/hyper.G-.link.zscore.txt > hyper.G-.link.zscore.all.optimized.links.txt
sed -i $'1 i\\\nprobe\tgeneSymbol\tgeneID\tZ.real\tpval' hyper.G-.link.zscore.all.optimized.links.txt
cd ..
fi
cd ..
echo "step4 done"
fi
# step 5 #
if [[ $step5 == "T" ]];
then
mkdir step5
cd ./step5
echo "step5 started"
findhypoGpos=$(head -n 49 ../settings.txt | tail -n 1)
eval $findhypoGpos
if [[ $findhypoGpos == "T" ]];
then
Rscript ../scripts/step5/make.summary.of.results.real.hypometh.G+.source.FIN.R
fi
findhypoGneg=$(head -n 50 ../settings.txt | tail -n 1)
eval $findhypoGneg
if [[ $findhypoGneg == "T" ]];
then
Rscript ../scripts/step5/make.summary.of.results.real.hypometh.G-.source.FIN.R
fi
findhyperGpos=$(head -n 56 ../settings.txt | tail -n 1)
eval $findhyperGpos
if [[ $findhyperGpos == "T" ]];
then
Rscript ../scripts/step5/make.summary.of.results.real.hypermeth.G+.source.FIN.R
fi
findhyperGneg=$(head -n 57 ../settings.txt | tail -n 1)
eval $findhyperGneg
if [[ $findhyperGneg == "T" ]];
then
Rscript ../scripts/step5/make.summary.of.results.real.hypermeth.G-.source.FIN.R
fi
hypoGposHistogram=$(head -n 67 ../settings.txt | tail -n 1)
eval $hypoGposHistogram
if [[ $hypoGposHistogram == "T" ]];
then
Rscript ../scripts/step5/make.histogram.of.results.real.hypometh.G+.source.R
fi
hypoGnegHistogram=$(head -n 68 ../settings.txt | tail -n 1)
eval $hypoGnegHistogram
if [[ $hypoGnegHistogram == "T" ]];
then
Rscript ../scripts/step5/make.histogram.of.results.real.hypometh.G-.source.R
fi
hyperGposHistogram=$(head -n 69 ../settings.txt | tail -n 1)
eval $hyperGposHistogram
if [[ $hyperGposHistogram == "T" ]];
then
Rscript ../scripts/step5/make.histogram.of.results.real.hypermeth.G+.source.R
fi
hyperGnegHistogram=$(head -n 70 ../settings.txt | tail -n 1)
eval $hyperGnegHistogram
if [[ $hyperGnegHistogram == "T" ]];
then
Rscript ../scripts/step5/make.histogram.of.results.real.hypermeth.G-.source.R
fi
hypoGposScatter=$(head -n 74 ../settings.txt | tail -n 1)
eval $hypoGposScatter
if [[ $hypoGposScatter == "T" ]];
then 
Rscript ../scripts/step5/make.simple.scatterplots.hypometh.G+.source.R
fi
hypoGnegScatter=$(head -n 75 ../settings.txt | tail -n 1)
eval $hypoGnegScatter
if [[ $hypoGnegScatter == "T" ]];
then
Rscript ../scripts/step5/make.simple.scatterplots.hypometh.G-.source.R
fi
hyperGposScatter=$(head -n 76 ../settings.txt | tail -n 1)
eval $hyperGposScatter
if [[ $hyperGposScatter == "T" ]];
then
Rscript ../scripts/step5/make.simple.scatterplots.hypermeth.G+.source.R
fi
hyperGnegScatter=$(head -n 77 ../settings.txt | tail -n 1)
eval $hyperGnegScatter
if [[ $hyperGnegScatter == "T" ]];
then
Rscript ../scripts/step5/make.simple.scatterplots.hypermeth.G-.source.R
fi
hypoGposTracks=$(head -n 82 ../settings.txt | tail -n 1)
eval $hypoGposTracks
if [[ $hypoGposTracks == "T" ]];
then
Rscript ../scripts/step5/make.ucsc.bed.file.for.hypometh.G+.source.R
fi
hypoGnegTracks=$(head -n 83 ../settings.txt | tail -n 1)
eval $hypoGnegTracks
if [[ $hypoGnegTracks == "T" ]];
then
Rscript ../scripts/step5/make.ucsc.bed.file.for.hypometh.G-.source.R
fi
hyperGposTracks=$(head -n 84 ../settings.txt | tail -n 1)
eval $hyperGposTracks
if [[ $hyperGposTracks == "T" ]];
then
Rscript ../scripts/step5/make.ucsc.bed.file.for.hypermeth.G+.source.R
fi
hyperGnegTracks=$(head -n 85 ../settings.txt | tail -n 1)
eval $hyperGnegTracks
if [[ $hyperGnegTracks == "T" ]];
then
Rscript ../scripts/step5/make.ucsc.bed.file.for.hypermeth.G-.source.R
fi
hypoGposStates=$(head -n 90 ../settings.txt | tail -n 1)
eval $hypoGposStates
if [[ $hypoGposStates == "T" ]];
then
Rscript ../scripts/step5/find.tumor.states.for.links.hypometh.G+.source.realfast.R
cd ./hypo.G+.output.states
allstates=$(wc -l *links.states.table.txt)
allstatesc=$(echo $allstates | awk '{print $1}')
if [[ $allstatesc > 0 ]];
rm cg*states.txt;
then echo "hypoGpos states of samples are annotated";
fi
cd ..
fi
hypoGnegStates=$(head -n 91 ../settings.txt | tail -n 1)
eval $hypoGnegStates
if [[ $hypoGnegStates == "T" ]];
then
Rscript ../scripts/step5/find.tumor.states.for.links.hypometh.G-.source.realfast.R
cd ./hypo.G-.output.states
allstates=$(wc -l *links.states.table.txt)
allstatesc=$(echo $allstates | awk '{print $1}')
if [[ $allstatesc > 0 ]];
rm cg*states.txt;
then echo "hypoGneg states of samples are annotated";
fi
cd ..
fi
hyperGposStates=$(head -n 92 ../settings.txt | tail -n 1)
eval $hyperGposStates
if [[ $hyperGposStates == "T" ]];
then
Rscript ../scripts/step5/find.tumor.states.for.links.hypermeth.G+.source.realfast.R
cd ./hyper.G+.output.states
allstates=$(wc -l *links.states.table.txt)
allstatesc=$(echo $allstates | awk '{print $1}')
if [[ $allstatesc > 0 ]];
rm cg*states.txt;
then echo "hyperGpos states of samples are annotated";
fi
cd ..
fi
hyperGnegStates=$(head -n 93 ../settings.txt | tail -n 1)
eval $hyperGnegStates
if [[ $hyperGnegStates == "T" ]];
then
Rscript ../scripts/step5/find.tumor.states.for.links.hypermeth.G-.source.realfast.R
cd ./hyper.G-.output.states
allstates=$(wc -l *links.states.table.txt)
allstatesc=$(echo $allstates | awk '{print $1}')
if [[ $allstatesc > 0 ]];
rm cg*states.txt;
then echo "hyperGneg states of samples are annotated";
fi
cd ..
fi
makeScatter4probe=$(head -n 102 ../settings.txt | tail -n 1)
eval $makeScatter4probe
if [[ $makeScatter4probe == "T" ]];
then
Rscript ../scripts/step5/make.simple.scatterplots.select.probe.R
fi
# Survival Analysis #
hypoGposSurvival=$(head -n 105 ../settings.txt | tail -n 1)
eval $hypoGposSurvival
if [[ $hypoGposSurvival == "T" ]];
then
mkdir "hypo.G+.output.survival"
Rscript ../scripts/step5/make.survival.plots.hypometh.G+.R
fi
hypoGnegSurvival=$(head -n 106 ../settings.txt | tail -n 1)
eval $hypoGnegSurvival
if [[ $hypoGnegSurvival == "T" ]];
then
mkdir "hypo.G-.output.survival"
Rscript ../scripts/step5/make.survival.plots.hypometh.G-.R
fi
hyperGposSurvival=$(head -n 107 ../settings.txt | tail -n 1)
eval $hyperGposSurvival
if [[ $hyperGposSurvival == "T" ]];
then
mkdir "hyper.G+.output.survival"
Rscript ../scripts/step5/make.survival.plots.hypermeth.G+.R
fi
hyperGnegSurvival=$(head -n 108 ../settings.txt | tail -n 1)
eval $hyperGnegSurvival
if [[ $hyperGnegSurvival == "T" ]];
then
mkdir "hyper.G-.output.survival"
Rscript ../scripts/step5/make.survival.plots.hypermeth.G-.R
fi
# circos plots #
hypoGposCircos=$(head -n 117 ../settings.txt | tail -n 1)
eval $hypoGposCircos
if [[ $hypoGposCircos == "T" ]];
then
mkdir "hypo.G+.output.circos"
Rscript ../scripts/step5/make.circos.plots.hypometh.G+.R
fi
hypoGnegCircos=$(head -n 118 ../settings.txt | tail -n 1)
eval $hypoGnegCircos
if [[ $hypoGnegCircos == "T" ]];
then
mkdir "hypo.G-.output.circos"
Rscript ../scripts/step5/make.circos.plots.hypometh.G-.R
fi
hyperGposCircos=$(head -n 119 ../settings.txt | tail -n 1)
eval $hyperGposCircos
if [[ $hyperGposCircos == "T" ]];
then
mkdir "hyper.G+.output.circos"
Rscript ../scripts/step5/make.circos.plots.hypermeth.G+.R
fi
hyperGnegCircos=$(head -n 120 ../settings.txt | tail -n 1)
eval $hyperGnegCircos
if [[ $hyperGnegCircos == "T" ]];
then
mkdir "hyper.G-.output.circos"
Rscript ../scripts/step5/make.circos.plots.hypermeth.G-.R
fi
# TAD analysis #
hypoGposTAD=$(head -n 125 ../settings.txt | tail -n 1)
eval $hypoGposTAD
if [[ $hypoGposTAD == "T" ]];
then
mkdir "hypo.G+.output.TAD"
Rscript ../scripts/step5/make.TAD.tables.hypometh.G+.R
fi
hypoGnegTAD=$(head -n 126 ../settings.txt | tail -n 1)
eval $hypoGnegTAD
if [[ $hypoGnegTAD == "T" ]];
then
mkdir "hypo.G-.output.TAD"
Rscript ../scripts/step5/make.TAD.tables.hypometh.G-.R
fi
hyperGposTAD=$(head -n 127 ../settings.txt | tail -n 1)
eval $hyperGposTAD
if [[ $hyperGposTAD == "T" ]];
then
mkdir "hyper.G+.output.TAD"
Rscript ../scripts/step5/make.TAD.tables.hypermeth.G+.R
fi
hyperGnegTAD=$(head -n 128 ../settings.txt | tail -n 1)
eval $hyperGnegTAD
if [[ $hyperGnegTAD == "T" ]];
then
mkdir "hyper.G-.output.TAD"
Rscript ../scripts/step5/make.TAD.tables.hypermeth.G-.R
fi
# Methylation heat map analysis #
hypoGposMetHeatmap=$(head -n 133 ../settings.txt | tail -n 1)
eval $hypoGposMetHeatmap
if [[ $hypoGposMetHeatmap == "T" ]];
then
mkdir "hypo.G+.output.probe.heatmap"
Rscript ../scripts/step5/make.met.heatmap.hypometh.G+.R
fi
hypoGnegMetHeatmap=$(head -n 134 ../settings.txt | tail -n 1)
eval $hypoGnegMetHeatmap
if [[ $hypoGnegMetHeatmap == "T" ]];
then
mkdir "hypo.G-.output.probe.heatmap"
Rscript ../scripts/step5/make.met.heatmap.hypometh.G-.R
fi
hyperGposMetHeatmap=$(head -n 135 ../settings.txt | tail -n 1)
eval $hyperGposMetHeatmap
if [[ $hyperGposMetHeatmap == "T" ]];
then
mkdir "hyper.G+.output.probe.heatmap"
Rscript ../scripts/step5/make.met.heatmap.hypermeth.G+.R
fi
hyperGnegMetHeatmap=$(head -n 136 ../settings.txt | tail -n 1)
eval $hyperGnegMetHeatmap
if [[ $hyperGnegMetHeatmap == "T" ]];
then
mkdir "hyper.G-.output.probe.heatmap"
Rscript ../scripts/step5/make.met.heatmap.hypermeth.G-.R
fi
# complex scatterplots #
purityinfo=$(head -n 25 ../settings.txt | tail -n 1)
eval $purityinfo
SMdataset=$(head -n 39 ../settings.txt | tail -n 1)
eval $SMdataset
CNVdataset=$(head -n 40 ../settings.txt | tail -n 1)
eval $CNVdataset
if [[ $purityinfo == "T" && $SMdataset == "T" && $CNVdataset == "T" ]];
then
hypoGposCScatter=$(head -n 96 ../settings.txt | tail -n 1)
eval $hypoGposCScatter
if [[ $hypoGposCScatter == "T" ]];
then
Rscript ../scripts/step5/make.complex.scatterplots.hypometh.G+.source.R
fi
hypoGnegCScatter=$(head -n 97 ../settings.txt | tail -n 1)
eval $hypoGnegCScatter
if [[ $hypoGnegCScatter == "T" ]];
then
Rscript ../scripts/step5/make.complex.scatterplots.hypometh.G-.source.R
fi
hyperGposCScatter=$(head -n 98 ../settings.txt | tail -n 1)
eval $hyperGposCScatter
if [[ $hyperGposCScatter == "T" ]];
then
Rscript ../scripts/step5/make.complex.scatterplots.hypermeth.G+.source.R
fi
hyperGnegCScatter=$(head -n 99 ../settings.txt | tail -n 1)
eval $hyperGnegCScatter
if [[ $hyperGnegCScatter == "T" ]];
then
Rscript ../scripts/step5/make.complex.scatterplots.hypermeth.G-.source.R
fi
fi
if [[ $purityinfo == "T" && $SMdataset == "F" && $CNVdataset == "F" ]];
then
hypoGposCScatter=$(head -n 96 ../settings.txt | tail -n 1)
eval $hypoGposCScatter
if [[ $hypoGposCScatter == "T" ]];
then
Rscript ../scripts/step5/make.purity.scatterplots.hypometh.G+.source.R
fi
hypoGnegCScatter=$(head -n 97 ../settings.txt | tail -n 1)
eval $hypoGnegCScatter
if [[ $hypoGnegCScatter == "T" ]];
then
Rscript ../scripts/step5/make.purity.scatterplots.hypometh.G-.source.R
fi
hyperGposCScatter=$(head -n 98 ../settings.txt | tail -n 1)
eval $hyperGposCScatter
if [[ $hyperGposCScatter == "T" ]];
then
Rscript ../scripts/step5/make.purity.scatterplots.hypermeth.G+.source.R
fi
hyperGnegCScatter=$(head -n 99 ../settings.txt | tail -n 1)
eval $hyperGnegCScatter
if [[ $hyperGnegCScatter == "T" ]];
then
Rscript ../scripts/step5/make.purity.scatterplots.hypermeth.G-.source.R
fi
fi
if [[ $purityinfo == "F" && $SMdataset == "T" && $CNVdataset == "F" ]];
then
hypoGposCScatter=$(head -n 96 ../settings.txt | tail -n 1)
eval $hypoGposCScatter
if [[ $hypoGposCScatter == "T" ]];
then
Rscript ../scripts/step5/make.SM.scatterplots.hypometh.G+.source.R
fi
hypoGnegCScatter=$(head -n 97 ../settings.txt | tail -n 1)
eval $hypoGnegCScatter
if [[ $hypoGnegCScatter == "T" ]];
then
Rscript ../scripts/step5/make.SM.scatterplots.hypometh.G-.source.R
fi
hyperGposCScatter=$(head -n 98 ../settings.txt | tail -n 1)
eval $hyperGposCScatter
if [[ $hyperGposCScatter == "T" ]];
then
Rscript ../scripts/step5/make.SM.scatterplots.hypermeth.G+.source.R
fi
hyperGnegCScatter=$(head -n 99 ../settings.txt | tail -n 1)
eval $hyperGnegCScatter
if [[ $hyperGnegCScatter == "T" ]];
then
Rscript ../scripts/step5/make.SM.scatterplots.hypermeth.G-.source.R
fi
fi
if [[ $purityinfo == "F" && $SMdataset == "F" && $CNVdataset == "T" ]];
then
hypoGposCScatter=$(head -n 96 ../settings.txt | tail -n 1)
eval $hypoGposCScatter
if [[ $hypoGposCScatter == "T" ]];
then
Rscript ../scripts/step5/make.CNV.scatterplots.hypometh.G+.source.R
fi
hypoGnegCScatter=$(head -n 97 ../settings.txt | tail -n 1)
eval $hypoGnegCScatter
if [[ $hypoGnegCScatter == "T" ]];
then
Rscript ../scripts/step5/make.CNV.scatterplots.hypometh.G-.source.R
fi
hyperGposCScatter=$(head -n 98 ../settings.txt | tail -n 1)
eval $hyperGposCScatter
if [[ $hyperGposCScatter == "T" ]];
then
Rscript ../scripts/step5/make.CNV.scatterplots.hypermeth.G+.source.R
fi
hyperGnegCScatter=$(head -n 99 ../settings.txt | tail -n 1)
eval $hyperGnegCScatter
if [[ $hyperGnegCScatter == "T" ]];
then
Rscript ../scripts/step5/make.CNV.scatterplots.hypermeth.G-.source.R
fi
fi
if [[ $purityinfo == "F" && $SMdataset == "F" && $CNVdataset == "F" ]];
then 
hypoGposCScatter=$(head -n 96 ../settings.txt | tail -n 1)
eval $hypoGposCScatter
if [[ $hypoGposCScatter == "T" ]];
then
echo "please load purity estimates or somatic mutation or CNV datasets to make complex scatterplots"
echo "please make sure your parameters are correct in the settings"
fi
hypoGnegCScatter=$(head -n 97 ../settings.txt | tail -n 1)
eval $hypoGnegCScatter
if [[ $hypoGnegCScatter == "T" ]];
then
echo "please load purity estimates or somatic mutation or CNV datasets to make complex scatterplots"
echo "please make sure your parameters are correct in the settings"
fi
hyperGposCScatter=$(head -n 98 ../settings.txt | tail -n 1)
eval $hyperGposCScatter
if [[ $hyperGposCScatter == "T" ]];
then
echo "please load purity estimates or somatic mutation or CNV datasets to make complex scatterplots"
echo "please make sure your parameters are correct in the settings"
fi
hyperGnegCScatter=$(head -n 99 ../settings.txt | tail -n 1)
eval $hyperGnegCScatter
if [[ $hyperGnegCScatter == "T" ]];
then
echo "please load purity estimates or somatic mutation or CNV datasets to make complex scatterplots"
echo "please make sure your parameters are correct in the settings"
fi
fi
cd ..
echo "step5 done"
fi
echo "TENET processing complete."
