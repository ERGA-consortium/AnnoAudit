process OMARK {
    label 'process_single'
    label 'omark'
    containerOptions = '-B $PWD:$HOME'

    input:
    path(database)     
    path(omamer)
     
    output:
    path("omark_output"), emit: omark
    path("omark_output/_detailed_summary.txt"), emit: omark_results

    script:
    """
    # For Matplotlib:
    export MPLCONFIGDIR="\${PWD}/.matplotlib"
    mkdir -p "\${MPLCONFIGDIR}"

    # For Fontconfig:
    export FONTCONFIG_PATH="\${PWD}/.fontconfig"
    mkdir -p "\${FONTCONFIG_PATH}"

    # Set writable directory for ETE3
    export ETE3_CACHE_HOME="\${PWD}/.etetoolkit"
    mkdir -p "\${ETE3_CACHE_HOME}"

    mkdir omark_output
    omark -f ${omamer} -d ${database} -o omark_output
    """
}
