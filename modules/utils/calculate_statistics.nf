process CALCULATE_STATISTICS {
    label 'process_single'
    label 'biopython'

    input:
    path(ch_fasta)
    path(ch_gff)
    val(cds_only)
    
    output:
    path("statistics.tsv")             , emit: statistics
    path("intron_sequences.fasta")     , emit: intron_fasta
    path("isoforms_introns.bed")       , emit: intron_bed
    path("longest_isoform_introns.bed"), emit: longest_isoform_introns

    script:
    """
    python3 ${projectDir}/bin/annot_report.py ${ch_gff} ${ch_fasta}
    """
}
