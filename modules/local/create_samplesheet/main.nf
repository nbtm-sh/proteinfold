process CREATE_SAMPLESHEET_YAML {
    tag "$meta.id"
    label 'process_single'

    container 'docker://nbtmsh/samplesheet-utils'

    input:
    tuple val(meta), path(samplesheet)

    output:
    tuple val(meta), path("*.yaml"), emit: samplesheet

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    create-samplesheet \\
        --directory ./ \\
        --yaml \
        --output-file \$(sample-name --sanitise --index 0 ./*.fasta).yaml
    """

    stub:
    """
    echo "" > samplesheet.yaml
    """
}

