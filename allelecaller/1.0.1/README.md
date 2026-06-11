# Getting started

This repository includes the following processes : 

- AlleleCalling
- AlleleFiltering

## Use in nextflow

To use these processes in your pipeline, you must have built & pushed Docker image to an external registry.
Once the image available in your registry, you can use it in your nextflow process with the following scripts.

**AlleleCalling**

```groovy
process alleleCalling {

    publishDir "${params.output_dir}/${params.run_directory}/${(params.use_sample_id) ? id : hashed_id}/alleleCalling", mode: params.output_mode_in_use

    input:
        tuple val(hashed_id), val(id), val(genus), val(species), path(blast_kb), path(assembly)
        path(qc_kb), stageAs: "qc_kb"

    output:
        tuple val(hashed_id), path("outputs.json"), emit: outputs
        tuple val(hashed_id), path("stats_calls.json.gz"), emit: stats
        tuple val(hashed_id), path("allele_calls.xml.gz"), path("allele_calls.bam"), path("allele_calls.json.gz"), emit: allele_calls
        tuple val(hashed_id), path("calls_standard.json.gz"), emit: standard_calls
        tuple val(hashed_id), path("calls_core_standard.csv.gz"), path("calls_core_pcr.csv.gz"), emit: csv_core
        tuple val(hashed_id), path("calls_accessory_standard.csv.gz"), path("calls_accessory_pcr.csv.gz"), emit: csv_accessory


    script:
    """
    ngs-run \
    --sample-id $id \
    --publish-dir $task.publishDir.path \
    --assembly $assembly \
    --blast-kb.path $blast_kb \
    --qc-kb.path $qc_kb \
    --organism.genus $genus \
    ${species ? '--organism.species ' + species : ''} \
    --n-threads ${task.cpus}
    """

    stub:
    """
    ngs-run \
    --sample-id $id \
    --publish-dir $task.publishDir.path \
    --assembly $assembly \
    --blast-kb.path $blast_kb \
    --qc-kb.path $qc_kb \
    --organism.genus $genus \
    ${species ? '--organism.species ' + species : ''} \
    --n-threads ${task.cpus} \
    --stub
    """
}

``` 

**AlleleFiltering**
```groovy
process alleleFiltering {

    publishDir "${params.output_dir}/${params.run_directory}/${(params.use_sample_id) ? id : id}/alleleFiltering", mode: params.output_mode_in_use

    input:
        tuple val(hashed_id), val(id), val(genus), val(species), path(filtering_kb), path(assembly), path(alignment), path(calls_bam)
        path(qc_kb), stageAs: "qc_kb"

    output:
        tuple  val(hashed_id), path("outputs.json"), emit: outputs
        tuple  val(hashed_id), path("stats_calls.json.gz"), emit: stats
        tuple  val(hashed_id), path("allele_calls.bam"), path("allele_calls.json.gz"), emit: allele_calls
        tuple val(hashed_id), path("calls_standard.json.gz"), emit: standard_calls
        tuple  val(hashed_id), path("calls_core_standard.csv.gz"), path("calls_core_pcr.csv.gz"), emit: csv_core
        tuple  val(hashed_id), path("calls_accessory_standard.csv.gz"), path("calls_accessory_pcr.csv.gz"), emit: csv_accessory

    script:
    """
    ngs-run AlleleFiltering \
    --sample-id $id \
    --publish-dir $task.publishDir.path \
    --assembly $assembly \
    --alignment $alignment \
    --calls-bam $calls_bam \
    --filtering-kb.path $filtering_kb \
    --qc-kb.path $qc_kb \
    --organism.genus $genus \
    ${species ? '--organism.species ' + species : ''} \
    --n-threads ${task.cpus}
    """
    stub:
    """
    ngs-run AlleleFiltering \
    --sample-id $id \
    --publish-dir $task.publishDir.path \
    --assembly $assembly \
    --alignment $alignment \
    --calls-bam $calls_bam \
    --filtering-kb.path $filtering_kb \
    --qc-kb.path $qc_kb \
    --organism.genus $genus \
    ${species ? '--organism.species ' + species : ''} \
    --n-threads ${task.cpus} \
    --stub
    """
}
``` 

# Running locally
The following examples assume you have the required files locally.  
You can either call using the local source code:
```bash
poetry run ngs-run AlleleCalling --sample-id <sample_id> --publish-dir /tmp/ --assembly <assembly.fasta.gz> --blast-kb.path <blast_kb> --qc-kb.path <qc_kb> --organism.genus <organism.genus> --stub
poetry run ngs-run AlleleFiltering --sample-id <sample_id> --publish-dir /tmp/ --assembly <assembly.fasta.gz> --alignment <alignment.cram> --calls-bam <calls.bam> --filtering-kb.path <blast_kb> --qc-kb.path <qc_kb> --organism.genus <organism.genus> --stub
```

or using the docker image once it has been built:
```bash
docker run -v $(pwd):/app ngs-pipeline-process-allele-calling ngs-run  AlleleCalling --sample-id <sample_id> --publish-dir /tmp/ --assembly <assembly.fasta.gz> --blast-kb.path <blast_kb> --qc-kb.path <qc_kb>  --organism.genus <organism.genus>
docker run -v $(pwd):/app ngs-pipeline-process-allele-calling ngs-run  AlleleFiltering --sample-id <sample_id> --publish-dir /tmp/ --assembly <assembly.fasta.gz> --alignment <alignment.cram> --calls-bam <calls.bam> --filtering-kb.path <blast_kb> --qc-kb.path <qc_kb> --organism.genus <organism.genus>
```

When running locally you have to use the `--stub` flag as you are likely to be missing some dependencies.  
With the Docker image you can run both in stub and non-stub modes.