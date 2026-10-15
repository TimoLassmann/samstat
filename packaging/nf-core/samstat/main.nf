process SAMSTAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // FIXME after the bioconda release: the biocontainer tag (version--build hash), e.g. from
    // https://quay.io/repository/biocontainers/samstat?tab=tags
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samstat:2.4.0--FIXME' :
        'quay.io/biocontainers/samstat:2.4.0--FIXME' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.samstat.html")    , emit: html
    tuple val(meta), path("*.samstat.json")    , emit: json
    tuple val(meta), path("*.samstat_*_mqc.json"), optional: true, emit: mqc
    tuple val("${task.process}"), val('samstat'), eval("samstat -v | sed 's/SAMStat //'"), emit: versions_samstat, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // reports are named after the input file: stage it as <prefix>.<extension>
    def ext = reads.name.endsWith('.gz') ? reads.name.tokenize('.')[-2..-1].join('.') : reads.extension
    def input = "${prefix}.${ext}"
    """
    [ "${reads}" = "${input}" ] || ln -s "${reads}" "${input}"

    samstat \\
        $args \\
        --json \\
        -t $task.cpus \\
        ${input}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext = reads.name.endsWith('.gz') ? reads.name.tokenize('.')[-2..-1].join('.') : reads.extension
    def mqc = args.contains('--multiqc') ? "touch ${prefix}.${ext}.samstat_general_mqc.json" : ""
    """
    touch ${prefix}.${ext}.samstat.html
    echo '{}' > ${prefix}.${ext}.samstat.json
    ${mqc}
    """
}
