#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INFERCNV } from './modules/infercnv'

workflow {
    validateParams()

    def sample_id = params.sample_id ?: file(params.raw_counts_matrix).baseName
    def meta = [ id: sample_id ]

    Channel
        .of([
            meta,
            file(params.raw_counts_matrix, checkIfExists: true),
            file(params.annotations_file, checkIfExists: true),
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

workflow.onComplete {
    log.info ''
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Duration             : ${workflow.duration}"
    log.info "Success              : ${workflow.success}"
    log.info "Output directory     : ${params.outdir}"
}

def validateParams() {
    def required = [
        raw_counts_matrix: params.raw_counts_matrix,
        annotations_file : params.annotations_file,
        gene_order_file  : params.gene_order_file
    ]

    def missing = required.findAll { !it.value }.keySet()
    if (missing) {
        def formatted = missing.collect { "--${it}" }.join(', ')
        exit 1, "Missing required parameter(s): ${formatted}"
    }
}
