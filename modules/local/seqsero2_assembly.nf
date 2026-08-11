process SEQSERO2_ASSEMBLY {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqsero2=1.3.1"
    container "quay.io/biocontainers/seqsero2:1.3.1--pyhdfd78af_0"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_assembly_results.tsv"), emit: tsv
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    SeqSero2_package.py \\
        $args \\
        -p $task.cpus \\
        -d ${meta.id}_assembly_out \\
        -n ${meta.id} \\
        -i ${fasta}

    cp ${meta.id}_assembly_out/SeqSero_result.tsv ${meta.id}_assembly_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqsero2: \$(SeqSero2_package.py --version 2>&1 | sed 's/^.*SeqSero2 //; s/ .*\$//')
    END_VERSIONS
    """
}