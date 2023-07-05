# 1. Dog to others
# 2. Human to others
# 3. Mouse to others


## Dog to others --------------------------------------------
mkdir 01_split_state
mkdir 02_200bin
mkdir 03_liftover_min95_hg38
mkdir 03_liftover_min95_mm10
mkdir 04_ratio_mapped

# split state 
cd 00_input_dense
for bed in $(ls |grep bed);do
	pre=${bed%%.bed}
	echo $pre
	awk -v pre=$pre '{print $0 >> "../01_split_state/"pre"."$4".bed"}' $bed
done
cd ..

# bed to 200bp bins
cd 01_split_state
for bed in *.bed;do
	pre=${bed%%.bed}
	echo $pre
	bedtools makewindows -b $bed -w 200 > ../02_200bin/${pre}.200bp.bed
done
cd ..

# liftOver to human and mouse (-minMisMatch 0.95)
cd 02_200bin
count=1
for bed in *.bed;do
	pre=${bed%%.bed}
	echo $count $pre
	./liftOver $bed canFam3ToHg38.over.chain ../03_liftover_min95_hg38/${pre}.tohg38.bed ../03_liftover_min95_hg38/${pre}.tohg38.unmap.bed
	./liftOver $bed canFam3ToMm10.over.chain ../03_liftover_min95_mm10/${pre}.tomm10.bed ../03_liftover_min95_mm10/${pre}.tomm10.unmap.bed
	count=$(expr $count + 1)
done
cd ..

## Count and calculate ratio of mapped to hg38 peaks
cd 02_200bin
wc -l *.bed > ../04_ratio_mapped/01_Total_num_of_bin
cd ..

cd 03_liftover_min95_hg38
wc -l *.tohg38.bed > ../04_ratio_mapped/02_mapped_num_of_bin_hg38
cd ..

cd 03_liftover_min95_mm10
wc -l *.tomm10.bed > ../04_ratio_mapped/02_mapped_num_of_bin_mm10
cd ..


## Human to others --------------------------------------------
# 1. split to 200bp bins
# Selected tissue: E073 E106 E094 E086 E097 E066 E098 E096 E113
for bed in *.bed;do
	pre=${bed%%.*};echo $pre
	bedtools makewindows -b $bed -w 200 -i src > ${path}/${pre}.02kb.bed

# 2. liftover to canFam3 and mm10
for label in E073 E106 E094 E086 E097 E066 E098 E096 E113;do
	echo $label
	./liftOver 220215_Used/${label}_15_coreMarks_mnemonics.hg38.02kb.bed hg38ToCanFam3.over.chain.gz ../04_liftover_hg38_to_canFam3/${label}_15_coreMarks_mnemonics.hg38.02kb.toCf3.bed tmp1.unmap -bedPlus=4
	./liftOver 220215_Used/${label}_15_coreMarks_mnemonics.hg38.02kb.bed hg38ToMm10.over.chain.gz ../04_liftover_hg38_to_mm10/${label}_15_coreMarks_mnemonics.hg38.02kb.toMm10.bed tmp1.unmap -bedPlus=4
done

# 3. separate states to other files
for bed in *.bed;do
	pre=${bed%%_*}
	echo $pre;sed 's/\///g' $bed |awk -v prefix=$pre '{print $4 >> "separate_states/"prefix"_"$4}'
done

# 4. Counting
path='/media/vetbio/Extend_03/210301_Main_Ref_Epi/02_ChromHMM/220215_hg19_chromHMM_to_canFam3_mm10/05_Count'
cd /media/vetbio/Extend_03/210301_Main_Ref_Epi/02_ChromHMM/220215_hg19_chromHMM_to_canFam3_mm10/03_split_to_200bp_bins/220215_Used/separate_states
wc -l * |sed 's/_/ /' > ${path}/hg38_02kb_count.txt

cd /media/vetbio/Extend_03/210301_Main_Ref_Epi/02_ChromHMM/220215_hg19_chromHMM_to_canFam3_mm10/04_liftover_hg38_to_canFam3/separate_states
wc -l * |sed 's/_/ /' > ${path}/hg38_02kb_to_canFam3_map.txt

cd /media/vetbio/Extend_03/210301_Main_Ref_Epi/02_ChromHMM/220215_hg19_chromHMM_to_canFam3_mm10/04_liftover_hg38_to_mm10/separate_states
wc -l * |sed 's/_/ /' > ${path}/hg38_02kb_to_mm10_map.txt


## Mouse to others --------------------------------------------
# 1. split to 200bp bins
path="/media/vetbio/Extend_03/210301_Main_Ref_Epi/02_ChromHMM/220215_mm10_chromHMM_to_others_P0/03_split_to_200bp_bins"
for bed in *.bed;do
	pre=${bed%%.*};echo $pre
	bedtools makewindows -b $bed -w 200 -i src > ${path}/${pre}.02kb.bed
done

# 2. liftover to canFam3 and hg38
path="/media/vetbio/Extend_03/210301_Main_Ref_Epi/02_ChromHMM/220215_mm10_chromHMM_to_others_P0/04_liftover_mm10_to_canFam3"
for bed in *.bed;do
	pre=${bed%%.*};echo $pre
	./liftOver $bed mm10ToCanFam3.over.chain.gz ${path}/${pre}.mm10.02kb.cf3.bed ${path}/unmapped/${pre}.mm10.02kb.cf3.unmap
	./liftOver $bed mm10ToHg38.over.chain.gz ${path}/${pre}.mm10.02kb.hg38.bed ${path}/unmapped/${pre}.mm10.02kb.hg38.unmap.txt
done

# 3. separate states to other files
for bed in *.bed;do
	pre=${bed%%_*};echo $pre;sed 's/\///g' $bed |awk -v prefix=$pre '{print $4 >> "separate_states/"prefix"_"$4}'
done

## 4. Counting
path='/media/vetbio/Extend_03/210301_Main_Ref_Epi/02_ChromHMM/220215_mm10_chromHMM_to_others_P0/05_Count'
cd /media/vetbio/Extend_03/210301_Main_Ref_Epi/02_ChromHMM/220215_hg19_chromHMM_to_canFam3_mm10/03_split_to_200bp_bins/220215_Used/separate_states
wc -l * |sed 's/_/ /' > ${path}/mm10_02kb_count.txt

cd /media/vetbio/Extend_03/210301_Main_Ref_Epi/02_ChromHMM/220215_hg19_chromHMM_to_canFam3_mm10/04_liftover_hg38_to_canFam3/separate_states
wc -l * |sed 's/_/ /' > ${path}/hg38_02kb_to_canFam3_map.txt

cd /media/vetbio/Extend_03/210301_Main_Ref_Epi/02_ChromHMM/220215_hg19_chromHMM_to_canFam3_mm10/04_liftover_hg38_to_mm10/separate_states
wc -l * |sed 's/_/ /' > ${path}/hg38_02kb_to_mm10_map.txt

