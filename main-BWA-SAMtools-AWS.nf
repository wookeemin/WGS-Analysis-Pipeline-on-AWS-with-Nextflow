process {
  executor = 'awsbatch'
  queue = 'ngs-queue'              //생성한 AWS Batch Queue 이름
  container = 'biocontainers/bwa:v0.7.17_cv1'     // 사용할 Docker 이미지
  cpus = 4
  memory = '8 GB'
  time = '6h'

  // 로그나 중간 파일을 저장할 S3 경로
  scratch = true
  errorStrategy = 'retry'
}

aws {
  region = 'eu-central-1'          // 지역 설정 
  batch {
    cliPath = '/usr/local/bin/aws'             // CloudShell 또는 EC2에서 CLI 설치 위치
    jobRole = 'arn:aws:iam::123456789012:role/ecsTaskExecutionRole'
    volumes = 'my-ebs-vol'
    privileged = true
  }
  client {
    uploadStorageClass = 'STANDARD'
    maxRetries = 3
  }
}

workDir = 's3://my-ngs-bucket/work'       // Nextflow 작업 디렉토리
bucketDir = 's3://my-ngs-bucket/output'   // 결과 파일 저장 디렉토리

docker {
  enabled = true
}
