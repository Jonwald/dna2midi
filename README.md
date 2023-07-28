# RNA-SEQ / RNA-Exome pipleline
This pipeline has been built and configured for use in an AWS VPC with no open internet acess
 - Fastq files are pulled directly from basespace
 - Jobs are exectued in AWS batch
 - temporary files and results are output to S3 buckets
 - Containers are pulled from private ECR repositories
 - results are output back to basespace
