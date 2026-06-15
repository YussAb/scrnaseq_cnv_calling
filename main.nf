#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INFERCNV } from './modules/infercnv'
include { COPYKAT  } from './modules/copykat'

workflow {
    def run_infercnv = asBooleanParam(params.run_infercnv)
    def run_copykat = asBooleanParam(params.run_copykat)

    validateParams(run_infercnv, run_copykat)

    def sample_id = params.sample_id ?: file(params.raw_counts_matrix).baseName
    def meta = [ id: sample_id ]
    def raw_counts_matrix = file(params.raw_counts_matrix, checkIfExists: true)
    def annotations_file = file(params.annotations_file, checkIfExists: true)

    if (run_infercnv) {
        Channel
            .of([
                meta,
                raw_counts_matrix,
                annotations_file,
                file(params.gene_order_file, checkIfExists: true)
            ])
            .set { ch_infercnv_input }

        INFERCNV(
            ch_infercnv_input,
            params.ref_group_names,
            params.annotations_delim,
            params.cutoff,
            params.cluster_by_groups,
            params.denoise,
            params.hmm
        )
    }

    if (run_copykat) {
        Channel
            .of([ meta, raw_counts_matrix, annotations_file ])
            .set { ch_copykat_input }

        COPYKAT (
            ch_copykat_input,
            params.ref_group_names,
            params.annotations_delim,
            params.copykat_raw_counts_delim,
            params.copykat_id_type,
            params.copykat_cell_line,
            params.copykat_ngene_chr,
            params.copykat_low_dr,
            params.copykat_up_dr,
            params.copykat_win_size,
            params.copykat_ks_cut,
            params.copykat_distance,
            params.copykat_output_seg,
            params.copykat_plot_genes,
            params.copykat_genome,
            params.copykat_norm_cell_names
        )
    }
}

workflow.onComplete {
    log.info ''
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Duration             : ${workflow.duration}"
    log.info "Success              : ${workflow.success}"
    log.info "Output directory     : ${params.outdir}"
}

def validateParams(run_infercnv, run_copykat) {
    if (!run_infercnv && !run_copykat) {
        exit 1, 'At least one module must be enabled: --run_infercnv true or --run_copykat true'
    }

    def required = [
        raw_counts_matrix: params.raw_counts_matrix,
        annotations_file : params.annotations_file
    ]

    if (run_infercnv) {
        required.gene_order_file = params.gene_order_file
    }

    def missing = required.findAll { !it.value }.keySet()
    if (missing) {
        def formatted = missing.collect { "--${it}" }.join(', ')
        exit 1, "Missing required parameter(s): ${formatted}"
    }
}

def asBooleanParam(value) {
    if (value instanceof Boolean) {
        return value
    }

    if (value == null) {
        return false
    }

    value.toString().toLowerCase() in ['true', 't', '1', 'yes', 'y']
}
