/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowIsoseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.fasta, params.primers, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()


/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { PERL_BIOPERL }    from '../modules/local/perl/bioperl/main'    addParams( options: modules['PERL_BIOPERL']    )
include { GSTAMA_FILELIST } from '../modules/local/gstama/filelist/main' addParams( options: modules['GSTAMA_FILELIST'] )


/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
refine_params = ""
if (params.require_polya == true) { refine_params = "--require-polya --min-polya-length $params.min_polya_length" }

tama_collpase_params = ""
if      (params.capped == true)  { tama_collpase_params = "-x capped -b BAM -a ${params.five_prime} -m ${params.splice_junction} -z ${params.three_prime}" }
else if (params.capped == false) { tama_collpase_params = "-x no_cap -b BAM -a ${params.five_prime} -m ${params.splice_junction} -z ${params.three_prime}" }

pbccs_params      = " --min-rq ${params.rq}"
tama_merge_params = "-a ${params.five_prime} -m ${params.splice_junction} -z ${params.three_prime}"

include { PBCCS }                            from '../modules/nf-core/modules/pbccs/main'           addParams( options: [ args:"${pbccs_params}", publish_dir:"1_PBCCS" ]                                   )
include { LIMA }                             from '../modules/nf-core/modules/lima/main'            addParams( options: modules['LIMA']                                                                     )
include { ISOSEQ3_REFINE }                   from '../modules/nf-core/modules/isoseq3/refine/main'  addParams( options: [ args:"${refine_params}", publish_dir:"3_ISOSEQ3_REFINE" ]                         )
include { PBBAM_PBMERGE as PBBAM_PBMERGE_1 } from '../modules/nf-core/modules/pbbam/pbmerge/main'   addParams( options: modules['PBBAM_PBMERGE_1']                                                          )
include { ISOSEQ3_CLUSTER }                  from '../modules/nf-core/modules/isoseq3/cluster/main' addParams( options: modules['ISOSEQ3_CLUSTER']                                                          )
include { PBBAM_PBMERGE as PBBAM_PBMERGE_2 } from '../modules/nf-core/modules/pbbam/pbmerge/main'   addParams( options: modules['PBBAM_PBMERGE_2']                                                          )
include { SAMTOOLS_FASTQ }                   from '../modules/nf-core/modules/samtools/fastq/main'  addParams( options: modules['SAMTOOLS_FASTQ']                                                           )
include { MINIMAP2_ALIGN }                   from '../modules/nf-core/modules/minimap2/align/main'  addParams( options: modules['MINIMAP2_ALIGN']                                                           )
include { GUNZIP }                           from '../modules/nf-core/modules/gunzip/main'          addParams( options: modules['GUNZIP']                                                                   )
include { ULTRA_PIPELINE }                   from '../modules/nf-core/modules/ultra/pipeline/main'  addParams( options: modules['ULTRA']                                                                    )
include { SAMTOOLS_SORT }                    from '../modules/nf-core/modules/samtools/sort/main'   addParams( options: modules['SAMTOOLS_SORT']                                                            )
include { BAMTOOLS_SPLIT }                   from '../modules/nf-core/modules/bamtools/split/main'  addParams( options: modules['BAMTOOLS_SPLIT']                                                           )
include { GSTAMA_COLLAPSE }                  from '../modules/nf-core/modules/gstama/collapse/main' addParams( options: [ args:"${tama_collpase_params}", publish_dir:"11_GSTAMA_COLLAPSE", suffix: "_tc" ] )
include { GSTAMA_MERGE }                     from '../modules/nf-core/modules/gstama/merge/main'    addParams( options: [ args:"${tama_merge_params}",    publish_dir:"12_GSTAMA_MERGE" ]                   )
include { MULTIQC }                          from '../modules/nf-core/modules/multiqc/main'         addParams( options: multiqc_options                                                                     )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ISOSEQ {
    //
    // SET UP CHANNELS
    //

    ch_versions = Channel.empty()
    // Check if input directory contains bam and pbis files and set up channels
    if (params.input) {
        Channel // Check presence of bam files
        .fromPath(params.input + '/*.bam')
        .ifEmpty { exit 1, "Cannot find any bam(s) in input directory: ${params.input}\nNB: File names must finish by '.bam'\nNB: Path needs to be enclosed in quotes!" }
        Channel // Check presence of bam.pbi files
        .fromPath(params.input + '/*.bam.pbi')
        .ifEmpty { exit 1, "Cannot find any pbi(s) in input directory: ${params.input}\nNB: File names must finish by '.bam.pbi'\nNB: Path needs to be enclosed in quotes!" }


        Channel // Prepare channel for bams: duplicates item by num of chunk and extract id
        .fromPath(params.input + '/*.bam')
        .flatMap {  // Duplicate each samples to match to the number of chunks
            def array = []
            for ( i = 1 ; i <= params.chunk ; i++ ) { array << it }
            return array
        }
        .map {  // Set channel to [meta, bam, pbi]
            [
                [ id:it.toString().replaceAll(/.*\\//, '').replaceAll(/.bam$/, '') ],
                file(it),
                file(it.toString().replaceAll(/.bam$/, '.bam.pbi'))
            ]
        }
        .set { ch_pbccs_in }
    } else { exit 1, 'OPTION ERROR: bam/bam.pbi directory not provided or cannot be found!' }


    Channel // Prepare the pbccs chunk_num channel
        .from((1..params.chunk).step(1).toList())
        .set { ch_chunk_num }


    Channel // Prepare the pbccs chunk_num value channel
        .value(params.chunk)
        .set { ch_chunk_on }


    if (params.primers) {
        Channel // Prepare value channel for primers used for the library preparation
            .value(file(params.primers))
            .set { ch_primers }
    } else { exit 1, 'OPTION ERROR: primers file not provided or cannot be found!' }


    if (params.fasta) {
        Channel // Prepare value channel for reference genome fasta file => minimap2/uLTRA
            .value(file(params.fasta))
            .set { ch_fasta }
    } else { exit 1, 'OPTION ERROR: fasta file not provided or cannot be found' }


    if (params.ultra == true) {
        Channel // --> Prepare gtf value channel for ultra
            .value(file(params.gtf))
            .set { ch_gtf }
    } else { exit 1, "OPTION ERROR: gtf file not provided or cannot be found" }


    //
    // START PIPELINE
    //
    PBCCS(ch_pbccs_in, ch_chunk_num, params.chunk) // Generate CCS from raw reads

    PBCCS.out.bam // Update meta, update id (+chunkX) and store former id
    .map {
        def chk       = (it[1] =~ /.*\.(chunk\d+)\.bam/)[ 0 ][ 1 ]
        def id_former = it[0].id
        def id_new    = it[0].id + "." + chk
        return [ [id:id_new, id_former:id_former], it[1] ]
    }
    .set { ch_pbccs_bam_updated }

    LIMA(ch_pbccs_bam_updated, ch_primers)   // Remove primers from CCS
    ISOSEQ3_REFINE(LIMA.out.bam, ch_primers) // Discard CCS without polyA tails, remove it from the other

    ISOSEQ3_REFINE.out.bam                   // Group BAM (refined CCS) by sample (former id)
        .map { [ it[0].id_former , it[1] ] }
        .groupTuple(by: 0, size: params.chunk)
        .map { [ [ id:it[0] ], it[1] ] }
        .set { ch_grouped_bams }

    PBBAM_PBMERGE_1(ch_grouped_bams) // Merged BAM files into one  file

    // The user have the choice to run or not the cluster step
    if (params.run_cluster == true) {            // Run cluster step
        ISOSEQ3_CLUSTER(PBBAM_PBMERGE_1.out.bam) // Cluster CCS, output consensus of CCS ("transcripts") and not clustered CCS ("singletons")
        ISOSEQ3_CLUSTER.out.bam                  // Concat bam and singletons bam channels together and group bams using id
            .concat(ISOSEQ3_CLUSTER.out.singletons_bam)
            .groupTuple()
            .set { ch_cluster_and_singletons }

        PBBAM_PBMERGE_2(ch_cluster_and_singletons) // Merge transcripts and singletons bams files into one

        PBBAM_PBMERGE_2.out.bam // Add single_end option to meta
            .map { [ [ id:it[0].id, single_end:true ], it[1] ] }
            .set { ch_cluster_updated }
    }
    else { // do not run cluster step
        PBBAM_PBMERGE_1.out.bam  // Add single_end option to meta
            .map { [ [ id:it[0].id, single_end:true ], it[1] ] }
            .set { ch_cluster_updated }
    }

    SAMTOOLS_FASTQ(ch_cluster_updated) // Convert BAM into FASTQ

    // Align CCS (no cluster path) or singletons + transcripts (cluster path)
    // User can choose between minimap2 and uLTRA aligners
    if (params.ultra == true) {
        GUNZIP(SAMTOOLS_FASTQ.out.fastq)
        ULTRA_PIPELINE(GUNZIP.out.gunzip, ch_fasta, ch_gtf)
        PERL_BIOPERL(ULTRA_PIPELINE.out.sam) // Remove remove reads ending with GAP (N) in CIGAR string
    }
    else {
        MINIMAP2_ALIGN(SAMTOOLS_FASTQ.out.fastq, ch_fasta) // Align read against genome
        PERL_BIOPERL(MINIMAP2_ALIGN.out.paf)               // Remove remove reads ending with GAP (N) in CIGAR string
    }

    SAMTOOLS_SORT(PERL_BIOPERL.out.out)   // Sort and convert sam to bam
    BAMTOOLS_SPLIT(SAMTOOLS_SORT.out.bam) // Split aligned sequences by chromosomes

    BAMTOOLS_SPLIT.out.bam // Discard unmapped sequences and store the sample id into id_former
        .map { discard_unmapped_and_save_id(it) }
        .flatten()
        .buffer( size: 2 )
        .set { ch_tc_input }

    GSTAMA_COLLAPSE(ch_tc_input, ch_fasta) // Clean gene models

    GSTAMA_COLLAPSE.out.bed // replace id with the former sample id and group files by sample
    .map { [ [id:it[0].id_former], it[1] ] }
    .groupTuple()
    .set { ch_tcollapse }

    GSTAMA_FILELIST(ch_tcollapse, Channel.value("capped"), Channel.value("1,1,1")) // Generate the filelist file needed by TAMA merge

    ch_tcollapse // Synchronized bed files produced by TAMA collapse with file list file generated by GSTAMA_FILELIST
    .join( GSTAMA_FILELIST.out.tsv )
    .set { ch_tmerge_in }

    GSTAMA_MERGE(ch_tmerge_in.map { [ it[0], it[1] ] }, ch_tmerge_in.map { it[2] }) // Merge bed files from one sample into one

    //
    // MODULE: Pipeline reporting
    //
    ch_versions = ch_versions.mix(PBCCS.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(LIMA.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(ISOSEQ3_REFINE.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(PBBAM_PBMERGE_1.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(PBBAM_PBMERGE_2.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(ISOSEQ3_CLUSTER.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(BAMTOOLS_SPLIT.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(GSTAMA_COLLAPSE.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(GSTAMA_MERGE.out.versions.first().ifEmpty(null))

    if (params.ultra == true) {
        ch_versions = ch_versions.mix(ULTRA_PIPELINE.out.versions.first().ifEmpty(null))
    }
    else {
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first().ifEmpty(null))
    }

    ch_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_versions }

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowIsoseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
//     ch_multiqc_files = ch_multiqc_files.mix(PBCCS.out.report_json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(LIMA.out.summary.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(LIMA.out.counts.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
     )
     multiqc_report = MULTIQC.out.report.toList()
     ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}


/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}


/*
========================================================================================
    UTILS
========================================================================================
*/
def discard_unmapped_and_save_id(List row) {
    def array = []
    for ( i in row[1] ) {
        def seq = (i =~ /.*\.(REF_.+)\.bam/)[ 0 ][ 1 ]
        if (seq != "REF_unmapped") {
            def id_former = row[0].id
            def id_new    = row[0].id + "." + seq
            array <<  [ [id:id_new, id_former:id_former], i ]
        }
        else { } //println("FOUND UNMAPPED: " + i)
    }
    return array
}


/*
========================================================================================
    THE END
========================================================================================
*/
