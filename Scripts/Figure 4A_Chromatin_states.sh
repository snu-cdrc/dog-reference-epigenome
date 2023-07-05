in_fold='01_Chip_2narrow_3broad'
in_file=${in_fold}/'01_ChromHMM_configure_wNum.txt'
in_chromsize='canFam3.chrom.sizes'

out_fold_bin='02_binarized'
out_fold='03_Results'
out_fold_final='output'

in_javaMem='100G'
in_thread=24

#---------------------------------------
mkdir ${out_fold_bin}
mkdir ${out_fold}

# Binarize Bam
date '+%F  %r'
echo 'Binarize Bed'

java -jar -Xmx${in_javaMem} ChromHMM/ChromHMM.jar BinarizeBed -peaks -p ${in_thread} -b 200 ${in_chromsize} ${in_fold} ${in_file} ${out_fold_bin}

# Learn model
for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20;do
    date '+%F  %r'
    echo 'Learn model : '${i}
    echo "java -jar -Xmx${in_javaMem} ChromHMM/ChromHMM.jar LearnModel -holdcolumnorder -p ${in_thread} -b 200 ${out_fold_bin} ${out_fold}/${out_fold_final}_${i} ${i} canFam3"
    java -jar -Xmx${in_javaMem} ChromHMM/ChromHMM.jar LearnModel -holdcolumnorder -p ${in_thread} -b 200 ${out_fold_bin} ${out_fold}/${out_fold_final}_${i} ${i} canFam3
done
date '+%F  %r'


# Compare model
for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20;do 
	java -jar -Xmx${in_javaMem} ChromHMM/ChromHMM.jar CompareModels 03_Results/output_20/emissions_20.txt 03_Results/output_${i} 04_Compare/compare_ref20_w${i}
done


# Reorder states -------------------------------------
java -jar -Xmx${in_javaMem} ChromHMM/ChromHMM.jar Reorder -holdcolumnorder -o 05_Reordering/210504_input_reorder 03_Results/output_13/model_13.txt 05_Reordering/

# Make Segmant.bed file
java -jar -Xmx${in_javaMem} ChromHMM/ChromHMM.jar MakeSegmentation -b 200 05_Reordering/model_13.txt 02_binarized/ 05_Reordering/

# NeighborhoodEnrichment TSS and TES
for tissue in CL CR CO KI LI LU MG OV PA SP ST;do echo $tissue;java -jar ChromHMM/ChromHMM.jar NeighborhoodEnrichment 05_Reordering/${tissue}_13_segments.bed 06_Enrich_plot/annotation/01_Ensembl/Ensemble_UCSC_3_TSS_btmg.anchor 05_Reordering/${tissue}_13_neighborhoodEnrich_TSS;done
for tissue in CL CR CO KI LI LU MG OV PA SP ST;do echo $tissue;java -jar ChromHMM/ChromHMM.jar NeighborhoodEnrichment 05_Reordering/${tissue}_13_segments.bed 06_Enrich_plot/annotation/01_Ensembl/Ensemble_UCSC_3_TES_btmg.anchor 05_Reordering/${tissue}_13_neighborhoodEnrich_TES;done

# Fix color of states 
for tissue in CL CR CO KI LI LU MG OV PA SP ST;do
	echo $tissue
	awk -v OFS='\t' '{ if($4=="U1") print $1, $2, $3, "1\t0\t.", $2, $3, "234,34,34" }' 05_Reordering/${tissue}_13_segments.bed >> 05_Reordering/${tissue}_13_dense_tmp.bed
	awk -v OFS='\t' '{ if($4=="U2") print $1, $2, $3, "2\t0\t.", $2, $3, "255,102,0" }' 05_Reordering/${tissue}_13_segments.bed >> 05_Reordering/${tissue}_13_dense_tmp.bed
	awk -v OFS='\t' '{ if($4=="U3") print $1, $2, $3, "3\t0\t.", $2, $3, "255,192,0" }' 05_Reordering/${tissue}_13_segments.bed >> 05_Reordering/${tissue}_13_dense_tmp.bed
	awk -v OFS='\t' '{ if($4=="U4") print $1, $2, $3, "4\t0\t.", $2, $3, "212,186,56" }' 05_Reordering/${tissue}_13_segments.bed >> 05_Reordering/${tissue}_13_dense_tmp.bed
	awk -v OFS='\t' '{ if($4=="U5") print $1, $2, $3, "5\t0\t.", $2, $3, "0,112,192" }' 05_Reordering/${tissue}_13_segments.bed >> 05_Reordering/${tissue}_13_dense_tmp.bed
	awk -v OFS='\t' '{ if($4=="U6") print $1, $2, $3, "6\t0\t.", $2, $3, "85,142,213" }' 05_Reordering/${tissue}_13_segments.bed >> 05_Reordering/${tissue}_13_dense_tmp.bed
	awk -v OFS='\t' '{ if($4=="U7") print $1, $2, $3, "7\t0\t.", $2, $3, "142,180,227" }' 05_Reordering/${tissue}_13_segments.bed >> 05_Reordering/${tissue}_13_dense_tmp.bed
	awk -v OFS='\t' '{ if($4=="U8") print $1, $2, $3, "8\t0\t.", $2, $3, "119,147,60" }' 05_Reordering/${tissue}_13_segments.bed >> 05_Reordering/${tissue}_13_dense_tmp.bed
	awk -v OFS='\t' '{ if($4=="U9") print $1, $2, $3, "9\t0\t.", $2, $3, "148,138,84" }' 05_Reordering/${tissue}_13_segments.bed >> 05_Reordering/${tissue}_13_dense_tmp.bed
	awk -v OFS='\t' '{ if($4=="U10") print $1, $2, $3, "10\t0\t.", $2, $3, "191,191,191" }' 05_Reordering/${tissue}_13_segments.bed >> 05_Reordering/${tissue}_13_dense_tmp.bed
	awk -v OFS='\t' '{ if($4=="U11") print $1, $2, $3, "11\t0\t.", $2, $3, "156,115,178" }' 05_Reordering/${tissue}_13_segments.bed >> 05_Reordering/${tissue}_13_dense_tmp.bed
	awk -v OFS='\t' '{ if($4=="U12") print $1, $2, $3, "12\t0\t.", $2, $3, "112,48,160" }' 05_Reordering/${tissue}_13_segments.bed >> 05_Reordering/${tissue}_13_dense_tmp.bed
	awk -v OFS='\t' '{ if($4=="U13") print $1, $2, $3, "13\t0\t.", $2, $3, "255,255,255" }' 05_Reordering/${tissue}_13_segments.bed >> 05_Reordering/${tissue}_13_dense_tmp.bed
	sort -k1,1 -k2,2n  05_Reordering/${tissue}_13_dense_tmp.bed > 05_Reordering/${tissue}_13_dense.bed
	rm 05_Reordering/${tissue}_13_dense_tmp.bed
done

# Overlap enrichment -------------------------------------
# 1. Ensemble genes & Repeats
for tissue in CL CR CO KI LI LU MG OV PA SP ST;do echo $tissue;java -jar -Xmx${in_javaMem} ChromHMM/ChromHMM.jar OverlapEnrichment 05_Reordering/${tissue}_13_segments.bed 06_Enrich_plot/annotation/01_Ensembl 05_Reordering/${tissue}_13_overEnrich;done

# 2. Expressed / Unexpressed genes
for tissue in CL CR CO KI LI LU MG OV PA SP ST;do echo $tissue;java -jar -Xmx${in_javaMem} ChromHMM/ChromHMM.jar OverlapEnrichment 05_Reordering/${tissue}_13_segments.bed 06_Enrich_plot/annotation/02_Exp_FPKM_01/${tissue} 05_Reordering/${tissue}_13_overEnrich_exp;done

# 3. ZNF genes
for tissue in CL CR CO KI LI LU MG OV PA SP ST;do echo $tissue;java -jar -Xmx${in_javaMem} ChromHMM/ChromHMM.jar OverlapEnrichment 05_Reordering/${tissue}_13_segments.bed 06_Enrich_plot/annotation/06_ZNF/canFam31_En_v102_ZNF_gene.bed 06_Enrich_plot/annotation/06_ZNF/ZNF_${tissue}_13_overEnrich;done



