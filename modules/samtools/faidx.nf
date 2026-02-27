process FAIDX {
	label 'process_low'
	label 'minimap2'

	input:
	file(ch_genome)
	
	output:
	path("*.fai"), emit: fai

	script:
	"""
	samtools faidx ${ch_genome}
	"""
}
