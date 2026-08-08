process DOWNLOAD_NIRVANA {
    label 'dragen'
    container { (workflow.stubRun || task.executor == 'local') ? null : "${task.ext.dragen_container}" }
    publishDir "$params.outdir/", saveAs: { filename -> filename == "versions.yml" ? null : filename }, mode:'copy'

    output:
    path "nirvana_annotation_data",    emit: nirvana_annotation_data

    script:
    def exe_path = "${task.ext.dragen_path}"
    def dragen_params = [
        task.ext.dragen_args                          ?: "",
        params.extra_dragen_args                      ?: "",
        params.nirvana_assembly                       ? "-r ${params.nirvana_assembly} --versions-config /opt/edico/resources/annotation/all_annotations_${params.nirvana_assembly}.json" : ""
    ].join(' ').trim()

    """
    cat > credentials.json << EOF
    {
    "ApiKey": "${task.ext.dragen_user}",
    "ApiSecret": "${task.ext.dragen_password}"
    }
    EOF
    mkdir nirvana_annotation_data
    ${exe_path}/share/nirvana/DataManager download --credentials-file credentials.json ${dragen_params} \\
    -d nirvana_annotation_data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(${exe_path}/bin/dragen --version | tail -n 1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def exe_path = "${task.ext.dragen_path}"

    def dragen_params = [
        task.ext.dragen_args                          ?: "",
        params.extra_dragen_args                      ?: "",
        params.nirvana_assembly                       ? "-r ${params.nirvana_assembly} --versions-config /opt/edico/resources/annotation/all_annotations_${params.nirvana_assembly}.json" : ""
    ].join(' ').trim()

    """
    cat > credentials.json << EOF
    {
    "ApiKey": "${task.ext.dragen_user}",
    "ApiSecret": "${task.ext.dragen_password}"
    }
    EOF
    mkdir nirvana_annotation_data
    echo ${exe_path}/share/nirvana/DataManager download --credentials-file credentials.json ${dragen_params} \\
     -d nirvana_annotation_data > dragen_command.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo stub-run)
    END_VERSIONS
    """
}