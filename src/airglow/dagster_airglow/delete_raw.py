import dagster as dg
from dagster_ncsa import S3ResourceNCSA


@dg.asset(
    ins={"unzip_chunked_archive": dg.AssetIn()},
)
def delete_raw(
        context: dg.AssetExecutionContext,
        unzip_chunked_archive,
        s3: S3ResourceNCSA):
    context.log.info(f"Deleting raw files for {unzip_chunked_archive['site']} on {unzip_chunked_archive['observation_date']}")  # NOQA E501
    s3_client = s3.get_client()
    bucket = dg.EnvVar("DEST_BUCKET").get_value()

    # Delete the raw files from the bucket
    for raw_file in unzip_chunked_archive["raw_files"]:
        context.log.info(f"Deleting raw file: {raw_file}")
        s3_client.delete_object(
            Bucket=bucket,
            Key=raw_file
        )

    # Delete the cloud cover files
    for cloud_cover_file in unzip_chunked_archive['cloud_cover_files']:
        context.log.info(f"Deleting cloud cover file: {cloud_cover_file}")
        s3_client.delete_object(
            Bucket=bucket,
            Key=cloud_cover_file
        )

    # Delete the log file
    context.log.info(f"Deleting log file: {unzip_chunked_archive['instrument_log_file']}")
    s3_client.delete_object(
        Bucket=bucket,
        Key=unzip_chunked_archive['instrument_log_file']
    )
