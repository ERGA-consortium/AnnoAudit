process BUSCO_ASSEM {
    label 'process_medium'
    label 'busco'

    input:
    path(protein)
    val(busco_lineage)

    output:
    path ("assembled")                              , emit: assembled
    path ("assembled/short_summary.specific.*.json"), emit: results
    path ("assembled/short_summary.specific.*.txt"), emit: results_txt

    script:
    if (params.busco_database && busco_lineage) {
        """
        busco -m genome -i ${protein} -o assembled -c ${task.cpus} -l ${busco_lineage} --offline --download_path ${params.busco_database}/${params.odb_version}
        """
    } else if (busco_lineage) {
        """
        busco -m genome -i ${protein} -o assembled -c ${task.cpus} -l ${busco_lineage}
        """
    } else {
        """
        busco -m genome -i ${protein} -o assembled --auto-lineage -c ${task.cpus}
        """
    }   
}
