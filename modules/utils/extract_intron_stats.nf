process EXTRACT_INTRON_STATS {
    label 'process_single'
    label 'biopython'

    input:
    path(ch_fasta)
    val(genetic_code)
    path(general_stats)
    
    output:
    path("all_statistics.tsv"), emit: statistics
    path("intron_plot_data.pkl"), emit: plot_data
    
    script:
    """
    python3 ${projectDir}/bin/extract_intron_stats.py ${ch_fasta} ${genetic_code}
    cat statistics.tsv stop_codon_statistics.tsv >> all_statistics.tsv
    """
}
