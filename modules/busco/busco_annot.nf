process BUSCO_ANNOT {
    label 'process_medium'
    label 'busco'

    input:
    path(protein)
    val(busco_lineage)

    output:
    path ("annotated")                              , emit: annotated
    path ("annotated/short_summary.specific.*.json"), emit: results
    path ("annotated/short_summary.specific.*.txt"),  emit: results_txt

    script:
    if (params.busco_database && busco_lineage) {
        """
        busco -m prot -i ${protein} -o annotated -c ${task.cpus} -l ${busco_lineage} --offline --download_path ${params.busco_database}/${params.odb_version}
        """
    } else if (busco_lineage) {
        """
        busco -m prot -i ${protein} -o annotated -c ${task.cpus} -l ${busco_lineage}
        """
    } else {
        """
        busco -m prot -i ${protein} -o annotated --auto-lineage -c ${task.cpus}
        """
    }   
}
