#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SCRNASEQ_CNV_CALLING } from './workflows/scrnaseq_cnv_calling'

workflow {
    SCRNASEQ_CNV_CALLING()
}

workflow.onComplete {
    log.info ''
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Duration             : ${workflow.duration}"
    log.info "Success              : ${workflow.success}"
    log.info "Output directory     : ${params.outdir}"
}
