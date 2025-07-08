# WGS-Analysis-Pipeline-on-AWS-with-Nextflow
AWS + Nextflow-based Whole Genome Sequencing (WGS) pipeline project template and practical guide

Project Structure:
[image](https://github.com/user-attachments/assets/df4a3d12-c432-4e5c-a4c7-bf7aba8aa234)


Key Files Description

1. main.nf (Main Pipeline Definition)
   
params.reads = "./data/fastq/sample_{1,2}.fastq.gz"
params.ref = "./data/ref/genome.fa"

workflow {
    bwa_mem_ch = bwa_mem(params.reads, params.ref)
    dedup_bam_ch = mark_duplicates(bwa_mem_ch)
    recal_bam_ch = bqsr(dedup_bam_ch, params.ref, params.known_sites)
    variants_ch = haplotype_caller(recal_bam_ch, params.ref)
    filtered_vcf_ch = variant_filtering(variants_ch)
}


2. aws.config (AWS Batch Executor Settings)

process.executor = 'awsbatch'

aws {
  region = 'eu-central-1'
  batch {
    cliPath = '/home/ec2-user/miniconda3/bin/aws'
    jobRole = 'arn:aws:iam::123456789012:role/BatchJobRole'
    jobQueue = 'NGSJobQueue'
    jobDefinition = 'NGSJobDefinition'
  }
}


3.Practical Tips Summary

| Item                   | Description                                                                                           |
| ---------------------- | ----------------------------------------------------------------------------------------------------- |
| 🔹 Use S3              | Store all input/output files in `s3://your-bucket/...` for improved I/O and cost management           |
| 🔹 Use Spot Instances  | AWS Batch with Spot Instances can save up to 70–80% of compute cost                                   |
| 🔹 Co-locate Resources | Keep `S3`, `Batch`, and compute resources in the same region (e.g., `eu-central-1`)                   |
| 🔹 Control Parallelism | Use directives like `maxForks`, `cpus`, and `memory` in Nextflow for resource tuning                  |
| 🔹 Use Docker          | Utilize official Docker images (e.g., `broadinstitute/gatk`, `biocontainers/bwa`) for reproducibility |


4. Downloadable Template

You can download the complete Nextflow project structure as a ZIP file here:      e.g.,     nextflow_wgs_pipeline_template.zip 



