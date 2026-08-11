/*
 * Simulate amplicons from an assembled genome by using cutadapt to extract
 * the regions bounded by each forward/reverse primer pair (linked adapters).
 *
 * primer_pairs is a TSV with columns: name<TAB>fwd_seq<TAB>rev_seq
 * cutadapt's linked-adapter mode (-g FWD...REV_revcomp) keeps only the
 * sequence *between* matched primers when --discard-untrimmed is set,
 * which is what we treat as the in-silico amplicon.
 */
process CUTADAPT_SIMULATE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::cutadapt=4.9"
    container "quay.io/biocontainers/cutadapt:4.9--py39hff71179_0"

    input:
    tuple val(meta), path(fasta)
    path primer_pairs

    output:
    tuple val(meta), path("${meta.id}_simulated_amplicons.fasta.gz"), emit: amplicons
    path "${meta.id}_cutadapt.log",                                   emit: log
    path "versions.yml",                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Build a linked-adapter argument list from the primer pair TSV.
    # For each pair we generate: -g NAME=FWD...REVCOMP(REV)
    # so cutadapt anchors on the forward primer and requires the
    # reverse-complemented reverse primer downstream.
    ADAPTERS=""
    while IFS=\$'\\t' read -r name fwd rev; do
        [ -z "\$name" ] && continue
        case "\$name" in \\#*) continue ;; esac
        revc=\$(echo "\$rev" | rev | tr 'ACGTacgtRYKMBDHVrykmbdhv' 'TGCAtgcaYRMKVHDByrmkvhdb')
        ADAPTERS="\$ADAPTERS -g \${name}=\${fwd}...\${revc}"
    done < ${primer_pairs}

    cutadapt \\
        \$ADAPTERS \\
        --error-rate ${params.cutadapt_error_rate} \\
        --overlap ${params.cutadapt_min_overlap} \\
        --minimum-length ${params.cutadapt_min_len} \\
        --discard-untrimmed \\
        --fasta \\
        -o ${meta.id}_simulated_amplicons.fasta.gz \\
        ${fasta} > ${meta.id}_cutadapt.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}