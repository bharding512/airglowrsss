import dagster as dg
from dagster_ncsa import S3ResourceNCSA


class DeleteRawConfig(dg.Config):
    site: str
    observation_date: str
    raw_files: list[str]
    cloud_cover_files: list[str]
    instrument_log_file: str


@dg.asset(
    deps=["unzip_chunked_archive"],
)
def delete_raw(
        context: dg.AssetExecutionContext,
        s3: S3ResourceNCSA):
    # Convert the metadata produced by the unzip_chunked_archive asset to a
    # DeleteRawConfig object.
    upstream_metadata = context.instance.get_latest_materialization_event(
        dg.AssetKey("unzip_chunked_archive")).asset_materialization.metadata
    config = DeleteRawConfig(**upstream_metadata['delete_raw_config'].data)
    context.log.info(f"Delete Raw config: {config}")

    context.log.info(f"Deleting raw files for {config.site} on {config.observation_date}")  # NOQA E501
    s3_client = s3.get_client()
    bucket = dg.EnvVar("DEST_BUCKET").get_value()

    # Delete the raw files from the bucket
    for raw_file in config.raw_files:
        context.log.info(f"Deleting raw file: {raw_file}")
        s3_client.delete_object(
            Bucket=bucket,
            Key=raw_file
        )

    # Delete the cloud cover files
    for cloud_cover_file in config.cloud_cover_files:
        context.log.info(f"Deleting cloud cover file: {cloud_cover_file}")
        s3_client.delete_object(
            Bucket=bucket,
            Key=cloud_cover_file
        )

    # Delete the log file
    context.log.info(f"Deleting log file: {config.instrument_log_file}")
    s3_client.delete_object(
        Bucket=bucket,
        Key=config.instrument_log_file
    )
