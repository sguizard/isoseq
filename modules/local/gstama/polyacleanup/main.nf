// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GSTAMA_POLYACLEANUP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gs-tama=1.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gs-tama:1.0.3--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gs-tama:1.0.3--hdfd78af_0"

    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_tama.fa")                   , emit: fasta
    tuple val(meta), path("*_tama_polya_flnc_report.txt"), emit: report
    tuple val(meta), path("*_tama_tails.fa")             , emit: tails
    path "versions.yml"                                  , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if( "$fasta" == "${prefix}.fasta" | "$fasta" == "${prefix}.fa" ) error "Input and output names are the same, use the suffix option to disambiguate"
    """
    tama_flnc_polya_cleanup.py \\
        -f $fasta \\
        -p ${prefix} \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: 'No_version_available'
    END_VERSIONS
    """
}
