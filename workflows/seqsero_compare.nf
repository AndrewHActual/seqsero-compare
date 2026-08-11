include { SEROTYPE_SOURCES   } from '../subworkflows/local/serotype_sources'
include { COMPARE_SEROTYPES  } from '../modules/local/compare_serotypes'
include { MULTIQC            } from '../modules/local/multiqc'

workflow SEQSERO_COMPARE {

    // ---- Validate required params ----
    if (!params.input)        { error "Please provide a samplesheet with --input" }
    if (!params.primer_fasta) { error "Please provide the primer-pair TSV with --primer_fasta" }

    ch_primer_pairs = Channel.value( file(params.primer_fasta, checkIfExists: true) )

    // ---- Parse samplesheet ----
    // Expected columns:
    //   sample, amplicon_fastq_1, amplicon_fastq_2, assembly_fasta
    ch_input = Channel
        .fromPath( params.input, checkIfExists: true )
        .splitCsv( header: true )
        .map { row ->
            def meta = [ id: row.sample ]
            def reads = [ file(row.amplicon_fastq_1, checkIfExists: true),
                          file(row.amplicon_fastq_2, checkIfExists: true) ]
            def assembly = file(row.assembly_fasta, checkIfExists: true)
            return [ meta, reads, assembly ]
        }

    ch_amplicon_reads = ch_input.map { meta, reads, assembly -> [ meta, reads ] }
    ch_assembly_fasta = ch_input.map { meta, reads, assembly -> [ meta, assembly ] }

    // ---- Run all three serotyping sources ----
    SEROTYPE_SOURCES (
        ch_amplicon_reads,
        ch_assembly_fasta,
        ch_primer_pairs
    )

    // ---- Collect all result TSVs and compare ----
    COMPARE_SEROTYPES (
        SEROTYPE_SOURCES.out.amplicon.map  { meta, tsv -> tsv }.collect(),
        SEROTYPE_SOURCES.out.assembly.map  { meta, tsv -> tsv }.collect(),
        SEROTYPE_SOURCES.out.simulated.map { meta, tsv -> tsv }.collect()
    )

    // ---- Aggregate everything MultiQC-relevant ----
    ch_multiqc_files = Channel.empty()
        .mix( COMPARE_SEROTYPES.out.mqc )                               // custom-content table
        .mix( SEROTYPE_SOURCES.out.cutadapt_logs.collect().ifEmpty([]) ) // cutadapt module
        .collect()

    ch_multiqc_config = Channel.value( file(params.multiqc_config, checkIfExists: true) )

    MULTIQC (
        ch_multiqc_files,
        ch_multiqc_config,
        COMPARE_SEROTYPES.out.versions
    )
}