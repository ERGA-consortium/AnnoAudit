process BUSCO_PLOT {
    label 'process_low'
    label 'busco'

    input:
    path(busco_annot)
    path(busco_assem)

    output:
    path("busco_figure.png"), emit: busco_plot

    script:
    """
    busco --plot .
    """
}
