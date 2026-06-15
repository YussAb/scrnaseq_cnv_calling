process COPYCAT {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(raw_counts_matrix), path(annotations_file)
    val ref_group_names
    val annotations_delim
    val raw_counts_delim
    val id_type
    val cell_line
    val ngene_chr
    val low_dr
    val up_dr
    val win_size
    val ks_cut
    val distance
    val output_seg
    val plot_genes
    val genome
    val norm_cell_names

    output:
    tuple val(meta), path('copycat'), emit: results
    path 'versions.yml', emit: versions

    script:
    def ref_groups = ref_group_names instanceof List ? ref_group_names.join(',') : ref_group_names
    def normal_cells = norm_cell_names instanceof List ? norm_cell_names.join(',') : norm_cell_names
    def runner = "${projectDir}/modules/copycat/usr/bin/run_copycat.R"

    """
    Rscript "${runner}" \\
        --raw-counts-matrix "${raw_counts_matrix}" \\
        --annotations-file "${annotations_file}" \\
        --out-dir copycat \\
        --sample-id "${meta.id}" \\
        --ref-group-names "${ref_groups}" \\
        --norm-cell-names "${normal_cells}" \\
        --annotations-delim '${annotations_delim}' \\
        --raw-counts-delim '${raw_counts_delim}' \\
        --id-type "${id_type}" \\
        --cell-line "${cell_line}" \\
        --ngene-chr ${ngene_chr} \\
        --low-dr ${low_dr} \\
        --up-dr ${up_dr} \\
        --win-size ${win_size} \\
        --ks-cut ${ks_cut} \\
        --distance "${distance}" \\
        --output-seg ${output_seg} \\
        --plot-genes ${plot_genes} \\
        --genome "${genome}" \\
        --num-threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e "cat(as.character(getRversion()))")
        copykat: \$(Rscript -e "cat(as.character(utils::packageVersion('copykat')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p copycat
    touch copycat/${meta.id}.copykat.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: "stub"
        copykat: "stub"
    END_VERSIONS
    """
}
