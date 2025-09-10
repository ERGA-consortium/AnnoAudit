process GET_ANNOTATION {
    label 'process_single'
    label 'phylocontext'

    input:
    val(taxon_id)   

    output:
    path("results"), emit: phylocontext_out

    script:
    """
    get_annotations -t ${taxon_id} -r genus -o results
    """
}
