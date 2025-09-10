process PLOT_PSAURON {
    label 'process_single'
    label 'seaborn'

    input:
    path(ch_psauron_out)
    path(ch_predicted_proteome)
    
    output:
    path("PSAURON_plot.pdf"), emit: psauron_pdf
    path("PSAURON_plot.png"), emit: psauron_png

    script:
    """
    python3 ${projectDir}/bin/plot_psauron.py --input ${ch_psauron_out} --fasta ${ch_predicted_proteome}
    """
}
