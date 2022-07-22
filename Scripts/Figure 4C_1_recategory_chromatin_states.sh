## Merge canine dense bed
mkdir merge
mid='_13_dense';suff='200bp.bed';
for ti in CL CR CO KI LI LU MG OV PA SP ST;do echo $ti
	cat ${ti}${mid}.1.${suff} ${ti}${mid}.2.${suff} ${ti}${mid}.3.${suff} ${ti}${mid}.4.${suff} > merge/${ti}_cf3_pro.bed
	cat ${ti}${mid}.5.${suff} ${ti}${mid}.6.${suff} ${ti}${mid}.7.${suff} > merge/${ti}_cf3_enh.bed
	cat ${ti}${mid}.8.${suff} ${ti}${mid}.9.${suff} ${ti}${mid}.10.${suff} ${ti}${mid}.11.${suff} ${ti}${mid}.12.${suff} > merge/${ti}_cf3_het.bed
done

## Merge hg38 dense bed
mkdir merge
mid='_13_dense';suff='200bp.4col.tohg38.bed';
for ti in CL CR CO KI LI LU MG OV PA SP ST;do echo $ti
	cat ${ti}${mid}.1.${suff} ${ti}${mid}.2.${suff} ${ti}${mid}.3.${suff} ${ti}${mid}.4.${suff} > merge/${ti}_hg38_pro.bed
	cat ${ti}${mid}.5.${suff} ${ti}${mid}.6.${suff} ${ti}${mid}.7.${suff} > merge/${ti}_hg38_enh.bed
	cat ${ti}${mid}.8.${suff} ${ti}${mid}.9.${suff} ${ti}${mid}.10.${suff} ${ti}${mid}.11.${suff} ${ti}${mid}.12.${suff} > merge/${ti}_hg38_het.bed
done

## Merge mm10 dense bed
mkdir merge
mid='_13_dense';suff='200bp.4col.tomm10.bed';
for ti in CL CR CO KI LI LU MG OV PA SP ST;do echo $ti
	cat ${ti}${mid}.1.${suff} ${ti}${mid}.2.${suff} ${ti}${mid}.3.${suff} ${ti}${mid}.4.${suff} > merge/${ti}_mm10_pro.bed
	cat ${ti}${mid}.5.${suff} ${ti}${mid}.6.${suff} ${ti}${mid}.7.${suff} > merge/${ti}_mm10_enh.bed
	cat ${ti}${mid}.8.${suff} ${ti}${mid}.9.${suff} ${ti}${mid}.10.${suff} ${ti}${mid}.11.${suff} ${ti}${mid}.12.${suff} > merge/${ti}_mm10_het.bed
done
