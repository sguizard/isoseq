// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'


params.options = [:]
options        = initOptions(params.options)

process PERL_BIOPERL {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::perl-bioperl=1.7.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/perl-bioperl:1.7.2--pl526_11"
    } else {
        container "quay.io/biocontainers/perl-bioperl:1.7.2--pl526_11"
    }

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("${meta.id}${options.suffix}"), emit: out
    path "versions.yml"                                 , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    perl \\
        $options.args \\
        $file > $prefix

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( perl --version 2>&1|perl -ne 'print \$1 if (/(v\\d+\\.\\d+\\.\\d+)/)' )
    END_VERSIONS
    """
}
