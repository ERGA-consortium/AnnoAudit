process GET_TAXON_INFO {
    label 'process_single'
    label 'biopython'  

    input:
    val(email)
    val(taxon_id)
     
    output:
    path("Taxon_info.txt"), emit: taxon_info

    script:
    """
    python3 ${projectDir}/bin/get_taxon_info.py --email ${email} --taxon_id ${taxon_id}
    """
}
