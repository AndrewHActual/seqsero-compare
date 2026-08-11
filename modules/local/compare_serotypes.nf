process COMPARE_SEROTYPES {
    tag "serotype_comparison"
    label 'process_low'

    conda "conda-forge::pandas=2.2.2"
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    path amplicon_tsvs,  stageAs: 'amplicon/*'
    path assembly_tsvs,  stageAs: 'assembly/*'
    path simulated_tsvs, stageAs: 'simulated/*'

    output:
    path "serotype_comparison.tsv",     emit: table
    path "serotype_comparison.csv",     emit: csv
    path "seqsero2_comparison_mqc.tsv", emit: mqc
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    compare_serotypes.py \\
        --amplicon-dir  amplicon \\
        --assembly-dir  assembly \\
        --simulated-dir simulated \\
        --out-tsv serotype_comparison.tsv \\
        --out-csv serotype_comparison.csv \\
        --out-mqc seqsero2_comparison_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}