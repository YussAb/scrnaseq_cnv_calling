#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INFERCNV } from './modules/infercnv'
include { COPYCAT } from './modules/copycat'

workflow {
    def run_infercnv = asBooleanParam(params.run_infercnv)
    def run_copycat = asBooleanParam(params.run_copycat)

    validateParams(run_infercnv, run_copycat)

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

    if (run_copycat) {
        Channel
            .of([ meta, raw_counts_matrix, annotations_file ])
            .set { ch_copycat_input }

        COPYCAT(
            ch_copycat_input,
            params.ref_group_names,
            params.annotations_delim,
            params.copycat_raw_counts_delim,
            params.copycat_id_type,
            params.copycat_cell_line,
            params.copycat_ngene_chr,
            params.copycat_low_dr,
            params.copycat_up_dr,
            params.copycat_win_size,
            params.copycat_ks_cut,
            params.copycat_distance,
            params.copycat_output_seg,
            params.copycat_plot_genes,
            params.copycat_genome,
            params.copycat_norm_cell_names
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

def validateParams(run_infercnv, run_copycat) {
    if (!run_infercnv && !run_copycat) {
        exit 1, 'At least one module must be enabled: --run_infercnv true or --run_copycat true'
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
