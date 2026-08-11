process MULTIQC {
    label 'process_low'

    conda "bioconda::multiqc=1.25"
    container "quay.io/biocontainers/multiqc:1.25--pyhdfd78af_0"

    input:
    path multiqc_files, stageAs: "?/*"
    path multiqc_config
    path versions

    output:
    path "*multiqc_report.html", emit: report
    path "*_data",               emit: data
    path "*_plots",              emit: plots, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def config_arg  = multiqc_config ? "--config $multiqc_config" : ''
    """
    multiqc \\
        --force \\
        $config_arg \\
        $args \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed 's/^.*version //' )
    END_VERSIONS
    """
}