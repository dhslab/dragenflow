process ANNOTATE_SMALLVARIANTS {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-vep:release_105"
    publishDir "$params.outdir/${meta.id}/", saveAs: { filename -> filename.equals("versions.yml") ? null : filename }, mode:'copy'

    input:
    tuple val(meta), path(dragen_output)
    path reference
    path vepcache

    output:
    tuple val(meta), path("${meta.id}.hard-filtered.annotated.vcf.gz*", arity: '2'), emit: vcf
    path "versions.yml",    emit: versions

    script:
    """
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep \\
    --format vcf --vcf --fasta ${reference} --hgvs --symbol --term SO --flag_pick -o ${meta.id}.hard-filtered.annotated.vcf \\
    -i ${meta.id}.hard-filtered.vcf.gz --offline --cache --max_af --dir ${vepcache} && \\
    bgzip -c ${meta.id}.hard-filtered.annotated.vcf > ${meta.id}.hard-filtered.annotated.vcf.gz && \\
    tabix -p vcf ${meta.id}.hard-filtered.annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | cut -d ':' -f 2 | sed 's/\s*//g')
    END_VERSIONS
    /opt/vep/src/ensembl-vep/vep --dir ${vepcache} --show_cache_info | awk '{ print "    "\$1": "\$2; }' >> versions.yml
    """

    stub:
    """
    touch "${meta.id}.hard-filtered.annotated.vcf.gz"
    touch "${meta.id}.hard-filtered.annotated.vcf.gz.tbi"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    \$(cat $projectDir/assets/stub/versions/vep_version.yaml)
    END_VERSIONS
    """
}
