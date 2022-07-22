#!/usr/bin/

## Convert bed file to gff format -------------------------------------
for bed in *.bed;do 
	pre=$(echo $bed|cut -d '.' -f 1)
	echo $pre
	grep -v '#' $bed |awk '{OFS="\t"; print $1, $4, ".", $2, $3, ".",".",".", $4}' > ${pre}.gff
done

## SE calling (ROSE_main.py) ------------------------------------------
for tissue in cl cr co ki li lu mg ov pa sp st;do 
	chip=$(ls |grep -v '.bai' |grep '.bam' |grep 'k27ac'|grep $tissue'-')
	input=$(ls |grep -v '.bai' |grep '.bam'|grep 'input' |grep $tissue'-')
	peak=$(ls |grep '.gff' |grep ${tissue}'_')

	echo $chip" / "$input" / "$peak
	python2 ROSE_main.py \
	-g CANFAM3 \
	-i $peak \
	-r $chip \
	-c $input \
	-o 01_ROSE_output/$tissue \
	-s 12500 \
	-t 2000
done

## SE annotation (ROSE_geneMapper.py) ---------------------------------
for tissue in cl cr co ki li lu mg ov pa sp st;do
	input=$(ls |grep '_SuperEnhancers.table.txt' |grep $tissue'_')
	pre=$(echo $input |sed s/'_SuperEnhancers.table.txt'//)
	echo $pre
	python ROSE_geneMapper.py \
	-o $pre \
	-g CANFAM3 \
	-i $input
done


# Basic argument should be assigned with your directory and files
output_folder="/Directory/ROSE/output"
bed_file_name="Liver_K27ac"

## R restart (make pdf file)  -----------------------------------------
## outFolder /enhancerFile /enhancerName /wceName
R --no-save ${output_folder} ${output_folder}/${bed_file_name}_12KB_STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt ${bed_file_name} x < ROSE_callSuper.R

## Extract typical enhancer
bedtools intersect -v \
-a ${bed_file_name}.bed \
-b ${bed_file_name}_Gateway_SuperEnhancers.bed > ${bed_file_name}_Gateway_TypicalEnhancers.bed

## Extract constituent enhancer
bedtools intersect -wa \
-a ${bed_file_name}.bed \
-b ${bed_file_name}_Gateway_SuperEnhancers.bed > ${bed_file_name}_Gateway_Constituents.bed

## Total set annotation
python2 ROSE_geneMapper.py \
-o 191010_11tissue_SE \
-g CANFAM3 \
-i 191010_11tissue_SE.total.mg.4col_SuperEnhancers.table.txt
