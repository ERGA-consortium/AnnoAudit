process BEDTOOLS_COMPLEMENT {
	label 'process_low'
	label 'bedtools'

	input:
	path(annotation_gtf)
    path(genome_fai)

	output:
    path("intron_sorted.bed")  , emit: intron_bed

	script:
	"""
    cut -f1,2 ${genome_fai} > sizes.genome
    awk 'OFS="\t" {print \$1, \$2}' sizes.genome | sort -k1,1 > chromSizes.bed
    cat ${annotation_gtf} | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k4,4n -k5,5n"}' > in_sorted.gff
    bedtools complement -i in_sorted.gff -g chromSizes.bed 
    awk 'OFS="\t" {if (\$3 == "exon") print \$1, \$4-1, \$5}' in_sorted.gff > exon_sorted.bed
    bedtools complement -i <(cat exon_sorted.bed intergenic_sorted.bed | sort -k1,1 -k2,2n) -g chromSizes.bed > intron_sorted.bed
	"""
}
