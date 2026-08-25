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

    def dragen_hash_cmd = [
        params.build_rna_hash_table ? "--ht-build-rna-hashtable true" : "",
        params.build_hla_hash_table ? "--ht-build-hla-hashtable true" : "",
        params.build_methylation_hash_table ? "--ht-num-threads 40 --ht-seed-len 27 --ht-methylated-combined=true" : "",
        params.build_cnv_hash_table ? "--ht-build-cnv-hashtable true" : "",
        params.build_legacy_cnv_hash_table ? "--enable-cnv true" : ""
    ].findAll { it!= "" }.join(" ").trim()

    """
    mkdir dragen_hash
    ${exe_path}/bin/dragen --build-hash-table true --output-directory dragen_hash --ht-reference ${fasta} ${dragen_hash_cmd} ${args_license}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(${exe_path}/bin/dragen --version | tail -n 1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def exe_path = "${task.ext.dragen_path}"
    def args_license = (task.ext.dragen_user && task.ext.dragen_password) ? "--lic-server 'https://${task.ext.dragen_user}:${task.ext.dragen_password}@license.dragen.illumina.com'" : ""

    def dragen_hash_cmd = [
        params.build_rna_hash_table ? "--ht-build-rna-hashtable true" : "",
        params.build_hla_hash_table ? "--ht-build-hla-hashtable true" : "",
        params.build_methylation_hash_table ? "--ht-num-threads 40 --ht-seed-len 27 --ht-methylated-combined=true" : "",
        params.build_cnv_hash_table ? "--ht-build-cnv-hashtable true" : "",
        params.build_legacy_cnv_hash_table ? "--enable-cnv true" : ""
    ].findAll { it!= "" }.join(" ").trim()

    """
    mkdir dragen_hash
    echo ${exe_path}/bin/dragen --build-hash-table true --output-directory dragen_hash --ht-reference ${fasta} ${dragen_hash_cmd} ${args_license} > dragen_hash/dragen_command.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(echo stub-run)
    END_VERSIONS
    """
}