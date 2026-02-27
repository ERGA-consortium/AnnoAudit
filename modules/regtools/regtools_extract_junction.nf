process JUNCTION_EXTRACT {
	label 'process_low'
	label 'regtools'

	input:
	path(genome_bam)
	path(genome_bam_index)

	output:
    path("splice_junctions.bed")  , emit: junction_reads

	script:
	"""
	# Extract possible junction position from the BAM file
    regtools junctions extract -a 8 -m 50 -o splice_junctions.bed ${genome_bam} -s ${params.stranding}
	# Transform BED12 to BED6 format with corresponding intron position
	awk 'BEGIN{OFS="\\t"} 
         NF>=12 && \$6!="?" {
             split(\$11, bs, ","); 
             split(\$12, bst, ","); 
             if(bs[1]!="" && bs[2]!="") {
                 donor=\$2+bs[1];
                 acceptor=\$2+bst[2]-1;
                 print \$1, donor, acceptor+1, \$4, \$5, \$6;
             }
         }' splice_junctions.bed > regtools_true_introns.bed
	"""
}
