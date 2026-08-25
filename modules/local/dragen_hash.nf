process DRAGEN_HASH {
    label 'dragen'
    container { (workflow.stubRun || task.executor == 'local') ? null : "${task.ext.dragen_container}" }
    publishDir "$params.outdir/", saveAs: { filename -> filename == "versions.yml" ? null : filename }, mode:'copy'

    input:
    path (fasta)

    output:
    path "dragen_hash/*",    emit: dragen_hash

    script:
    def exe_path = "${task.ext.dragen_path}"
    def args_license = (task.ext.dragen_user && task.ext.dragen_password) ? "--lic-server 'https://${task.ext.dragen_user}:${task.ext.dragen_password}@license.dragen.illumina.com'" : ""

    """
    ${exe_path}/bin/dragen --build-hash-table true --output-directory dragen_hash --ht-reference ${fasta} --ht-build-rna-hashtable true --ht-build-hla-hashtable true --ht-methylated-cg true --ht-build-cnv-hashtable true ${args_license}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(${exe_path}/bin/dragen --version | tail -n 1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def exe_path = "${task.ext.dragen_path}"
    def args_license = (task.ext.dragen_user && task.ext.dragen_password) ? "--lic-server 'https://${task.ext.dragen_user}:${task.ext.dragen_password}@license.dragen.illumina.com'" : ""

    """
    mkdir dragen_hash
    echo ${exe_path}/bin/dragen --build-hash-table true --output-directory dragen_hash --ht-reference ${fasta} --ht-build-rna-hashtable true --ht-build-hla-hashtable true --ht-methylated-cg true --ht-build-cnv-hashtable true ${args_license} > dragen_hash/dragen_command.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo stub-run)
    END_VERSIONS
    """
}