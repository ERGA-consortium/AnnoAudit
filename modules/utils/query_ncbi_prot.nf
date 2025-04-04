process QUERY_NCBI_PROT {
    label 'process_single'
    label 'biopython'  

    input:
    val(email)
    val(taxon_id)
    val(target_count)
    val(batch_size)
     
    output:
    path("NCBI_protein_sequences.fasta"), emit: database
    path("species_counts.tsv")          , emit: species_count

    script:
    """
    python3 ${projectDir}/bin/query_protein_parallel.py --email ${email} --taxon_id ${taxon_id} --target_count ${target_count}
    """
}
