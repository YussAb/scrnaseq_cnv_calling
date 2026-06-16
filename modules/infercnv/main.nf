process INFERCNV {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(raw_counts_matrix), path(annotations_file), path(gene_order_file)
    val ref_group_names
    val annotations_delim
    val cutoff
    val cluster_by_groups
    val denoise
    val hmm
    val leiden_resolution
    val bayes_max_p_normal
    val plot_preliminary_cnv
    val plot_per_group

    output:
    tuple val(meta), path('infercnv'), emit: results
    path 'versions.yml', emit: versions

    script:
    def ref_groups = ref_group_names instanceof List ? ref_group_names.join(',') : ref_group_names
    def ref_group_args = ref_groups ? "--ref-group-names '${ref_groups}'" : ''
    //def runner = "${projectDir}/modules/infercnv/usr/bin/run_infercnv.R" Rscript "${runner}"

    """
    run_infercnv.R \\
        --raw-counts-matrix "${raw_counts_matrix}" \\
        --annotations-file "${annotations_file}" \\
        --gene-order-file "${gene_order_file}" \\
        --out-dir infercnv \\
        --annotations-delim '${annotations_delim}' \\
        ${ref_group_args} \\
        --cutoff ${cutoff} \\
        --cluster-by-groups ${cluster_by_groups} \\
        --denoise ${denoise} \\
        --hmm ${hmm} \\
        --leiden-resolution ${leiden_resolution} \\
        --bayes-max-p-normal ${bayes_max_p_normal} \\
        --plot-preliminary-cnv ${plot_preliminary_cnv} \\
        --plot-per-group ${plot_per_group} \\
        --num-threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e "cat(as.character(getRversion()))")
        infercnv: \$(Rscript -e "cat(as.character(utils::packageVersion('infercnv')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p infercnv
    touch infercnv/run.final.infercnv_obj

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: "stub"
        infercnv: "stub"
    END_VERSIONS
    """
}
