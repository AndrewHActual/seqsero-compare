#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SEQSERO_COMPARE } from './workflows/seqsero_compare'

workflow {
    SEQSERO_COMPARE()
}