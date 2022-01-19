process GSTAMA_FILELIST {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(bed)
    val cap
    val order

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for i in *.bed
    do
        echo -e "\${i}\\t${cap}\\t${order}\\t\${i}" >> ${prefix}.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gstama: \$( tama_merge.py -version | head -n1 )
    END_VERSIONS
    """
}
