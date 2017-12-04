#!/usr/bin/env bash

# Copy exome data to HDFS on a GCS cluster.

${GATK_HOME:-../..}/gatk-launch ParallelCopyGCSDirectoryIntoHDFSSpark \
    --input-gcs-path gs://broad-spark-eval-test-data/exome/ \
    --output-hdfs-directory hdfs://${GCS_CLUSTER}-m:8020/user/$USER/exome_spark_eval \
    -apiKey $API_KEY \
    -- \
    --sparkRunner GCS \
    --cluster $GCS_CLUSTER