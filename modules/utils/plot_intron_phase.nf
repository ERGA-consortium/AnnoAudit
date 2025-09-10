process PLOT_INTRON_PHASE {
    label 'process_single'
    label 'seaborn'

    input:
    path(plot_data)
    
    output:
    path('intron_distributions.pdf'), emit: intron_plot_pdf
    path('intron_distributions.png'), emit: intron_plot_png

    script:
    """
    python3 ${projectDir}/bin/plot_intron_phase.py --input_pkl ${plot_data}
    """
}
