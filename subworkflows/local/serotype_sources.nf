include { SEQSERO2_READS                        } from '../../modules/local/seqsero2_reads'
include { SEQSERO2_ASSEMBLY                     } from '../../modules/local/seqsero2_assembly'
include { CUTADAPT_SIMULATE                     } from '../../modules/local/cutadapt_simulate'
include { SEQSERO2_READS as SEQSERO2_SIMULATED  } from '../../modules/local/seqsero2_reads'

workflow SEROTYPE_SOURCES {
    take:
    ch_amplicon_reads   // [ meta, [reads] ]
    ch_assembly_fasta   // [ meta, fasta ]
    ch_primer_pairs     // path to primer-pair TSV

    main:
    ch_versions = Channel.empty()

    // Source 1: real amplicon reads
    SEQSERO2_READS ( ch_amplicon_reads )
    ch_versions = ch_versions.mix( SEQSERO2_READS.out.versions.first() )

    // Source 2: assembled genome fasta
    SEQSERO2_ASSEMBLY ( ch_assembly_fasta )
    ch_versions = ch_versions.mix( SEQSERO2_ASSEMBLY.out.versions.first() )

    // Source 3: in-silico amplicons cut from the assembly, then serotyped
    CUTADAPT_SIMULATE ( ch_assembly_fasta, ch_primer_pairs )
    ch_versions = ch_versions.mix( CUTADAPT_SIMULATE.out.versions.first() )

    SEQSERO2_SIMULATED ( CUTADAPT_SIMULATE.out.amplicons )
    ch_versions = ch_versions.mix( SEQSERO2_SIMULATED.out.versions.first() )

    emit:
    amplicon      = SEQSERO2_READS.out.tsv        // [ meta, tsv ]
    assembly      = SEQSERO2_ASSEMBLY.out.tsv     // [ meta, tsv ]
    simulated     = SEQSERO2_SIMULATED.out.tsv    // [ meta, tsv ]
    cutadapt_logs = CUTADAPT_SIMULATE.out.log     // path (for MultiQC)
    versions      = ch_versions
}