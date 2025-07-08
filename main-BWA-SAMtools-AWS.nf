process {
  executor = 'awsbatch'
  queue = 'ngs-queue'              //Created AWS Batch Queue name
  container = 'biocontainers/bwa:v0.7.17_cv1'     // Docker image
  cpus = 4
  memory = '8 GB'
  time = '6h'

  // S3 path to save log or intermediate files
  scratch = true
  errorStrategy = 'retry'
}

aws {
  region = 'eu-central-1'          // region
  batch {
    cliPath = '/usr/local/bin/aws'             // CLI installation path in CloudShell or EC2
  jobRole = 'arn:aws:iam::123456789012:role/ecsTaskExecutionRole'
    volumes = 'my-ebs-vol'
    privileged = true
  }
  client {
    uploadStorageClass = 'STANDARD'
    maxRetries = 3
  }
}

workDir = 's3://my-ngs-bucket/work'       // Nextflow work directory
bucketDir = 's3://my-ngs-bucket/output'   // output results saving directory

docker {
  enabled = true
}
