#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/isoseq
========================================================================================
 nf-core/isoseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/isoseq
----------------------------------------------------------------------------------------
*/

pbccs_params         = " --min-rq ${params.rq}"
lima_params          = " --isoseq --peek-guess"
refine_params        = ""
minimap2_params      = "-x splice:hq -uf --secondary=no -a"
tama_collpase_params = ""
tama_merge_params    = "-a ${params.five_prime} -m ${params.splice_junction} -z ${params.three_prime}"

if (params.require_polya == true) { refine_params = "--require-polya --min-polya-length $params.min_polya_length" }
if (params.capped == true)        { tama_collpase_params   = "-x capped -b BAM -a ${params.five_prime} -m ${params.splice_junction} -z ${params.three_prime}" }
else if (params.capped == false)  { tama_collpase_params   = "-x no_cap -b BAM -a ${params.five_prime} -m ${params.splice_junction} -z ${params.three_prime}" }

nextflow.enable.dsl=2
include { checkHostname }         from './lib/functions'
include { get_software_versions } from './lib/nf-core_process'
include { PBBAM_PBINDEX }         from './nf-core_modules/software/pbbam/pbindex/main.nf'
include { PBCCS }                 from './nf-core_modules/software/pbccs/main'           addParams( options: [ args:"${pbccs_params}" ] )
include { PBBAM_PBMERGE }         from './nf-core_modules/software/pbbam/pbmerge/main.nf'
include { LIMA }                  from './nf-core_modules/software/lima/main'            addParams( options: [ args:"${lima_params}" ] )
include { ISOSEQ3_REFINE }        from './nf-core_modules/software/isoseq3/refine/main'  addParams( options: [ args:" ${refine_params}" ] )
include { ISOSEQ3_CLUSTER }       from './nf-core_modules/software/isoseq3/cluster/main'
include { SAMTOOLS_FASTQ }        from './nf-core_modules/software/samtools/fastq/main'
include { MINIMAP2_ALIGN }        from './nf-core_modules/software/minimap2/align/main'  addParams( options: [ args:"${minimap2_params}" ] )
include { SAMTOOLS_SORTVIEW }     from './nf-core_modules/software/samtools/sortview/main'
include { BAMTOOLS_SPLIT }        from './nf-core_modules/software/bamtools/split/main'  addParams( options: [ args:"-reference" ] )
include { GSTAMA_COLLAPSE }       from './nf-core_modules/software/gstama/collapse/main' addParams( options: [ args:"${tama_collpase_params}" ] )
include { GSTAMA_MERGE }          from './nf-core_modules/software/gstama/merge/main' addParams( options: [ args:"${tama_merge_params}" ] )

log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/isoseq --input '*_R{1,2}.fastq.gz' -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}
////////////////////////////////////////////////////

/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

// Check if genome exists in the config file
// if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
//     exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(', ')}"
// }

if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) }

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}

// Stage config files
// ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
// ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

/*
 * Create a channel for input read files
 */

Channel
    .fromPath(params.input + '/*.bam')
    .map { row ->
        [
            row.toString().replaceAll(/.*\\//, '').replaceAll(/.bam$/, ''),
            file(row)
        ]
    }
    .ifEmpty { exit 1, "Cannot find any bam(s) in input directory: ${params.input}\nNB: Path needs to be enclosed in quotes!" }
    . set { ch_bams }

Channel
    .fromPath(params.input + '/*.bam.pbi')
    .map { row ->
        [
            row.toString().replaceAll(/.*\\//, '').replaceAll(/.bam.pbi$/, ''),
            file(row)
        ]
    }
    .ifEmpty { exit 1, "Cannot find any pbi(s) in input directory: ${params.input}\nNB: Path needs to be enclosed in quotes!" }
    . set { ch_pbis }

Channel
    .from((1..params.chunk).step(1).toList())
    .map { row -> [row, params.chunk] }
    .set { ch_chunk }

ch_bams
    .join(ch_pbis)
    .map { row ->
        [
            [id:row[0]],
            row[1],
            row[2]
        ]
    }
    .combine(ch_chunk)
    .set { ch_pbccs_in }

Channel
    .value(file(params.primers))
    .ifEmpty { exit 1, 'OPTION ERROR: primers file not provide or cannot be found' }
    .set { ch_primers }

Channel
    .value(file(params.fasta))
    .ifEmpty { exit 1, 'OPTION ERROR: fasta file not provide or cannot be found' }
    .set { ch_fasta }


////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = workflow.runName
// TODO nf-core: Report custom parameters here
summary['Input']            = params.input
summary['Fasta Ref']        = params.fasta
// summary['Data Type']        = params.single_end ? 'Single-End' : 'Paired-End'
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-isoseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/isoseq Workflow Summary'
    section_href: 'https://github.com/nf-core/isoseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
// process get_software_versions {
// }

/*
 * STEP 3 - Output Description HTML
 */
// process output_documentation {
//     publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode
//
//     input:
//     file output_docs
//     file images
//
//     output:
//     file 'results_description.html'
//
//     script:
//     """
//     markdown_to_html.py $output_docs -o results_description.html
//     """
// }


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/isoseq] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/isoseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/isoseq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
                mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/isoseq] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/isoseq]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/isoseq]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}


workflow {
    PBCCS(ch_pbccs_in)
    LIMA(PBCCS.out.bam, ch_primers)
    ISOSEQ3_REFINE(LIMA.out.bam, ch_primers)

    ISOSEQ3_REFINE.out.bam
        .map { row -> [ row[0].id , row[1] ] }
        .groupTuple(by: 0, size: params.chunk)
        .map { row -> [ [ id:row[0] ], row[1] ] }
        .set { ch_grouped_bams }

    PBBAM_PBMERGE(ch_grouped_bams)

    if (params.run_cluster == true) {
        ISOSEQ3_CLUSTER(PBBAM_PBMERGE.out.bam)

        ISOSEQ3_CLUSTER.out.bam
            .map { row -> [ [ id:row[0].id, single_end:true ], row[1] ] }
            .set { ch_cluster_updated }
    }
    else {
        PBBAM_PBMERGE.out.bam
            .map { row -> [ [ id:row[0].id, single_end:true ], row[1] ] }
            .set { ch_cluster_updated }
    }

    SAMTOOLS_FASTQ(ch_cluster_updated)
    MINIMAP2_ALIGN(SAMTOOLS_FASTQ.out.fastq, ch_fasta)
    SAMTOOLS_SORTVIEW(MINIMAP2_ALIGN.out.paf)
    BAMTOOLS_SPLIT(SAMTOOLS_SORTVIEW.out.bam)
    
    BAMTOOLS_SPLIT.out.bam
        .map { 
            out = []
            for ( i in it[1] ) {
                if (!(i =~ /unmapped/)) {out << [ it[0], i ]}
                else {
                    //println("FOUND UNMAPPED: " + i)
                }
            }
            return out
        }
        .flatten()
        .buffer( size: 2 )
        .set { ch_tc_input }
    
    //ch_tc_input.view()

    GSTAMA_COLLAPSE(ch_tc_input, ch_fasta)
    GSTAMA_MERGE(GSTAMA_COLLAPSE.out.bed.groupTuple())
    //GSTAMA_MERGE.out.bed.view()
}
