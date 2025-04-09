import dagster as dg
from botocore.exceptions import ClientError
from dagster import EnvVar, RunConfig
from dagster_ncsa import S3ResourceNCSA

from dagster_airglow.assets import ChunkedArchiveConfig


def list_files(bucket: str, prefix: str, s3_client) -> list[str]:
    if not prefix.endswith("/"):
        prefix += "/"

    try:
        # List all objects in the directory
        paginator = s3_client.get_paginator("list_objects_v2")
        objects = []

        # Use pagination to handle large directories
        for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
            if "Contents" in page:
                for obj in page["Contents"]:
                    objects.append(obj["Key"])
        return objects

    except ClientError as e:
        raise e


def group_files_by_date(file_list):
    date_files_dict = {}

    for file in file_list:
        # Extract the date part from the filename
        # The date appears in format YYYYMMDD
        filename = file.split('/')[-1]  # Get just the filename without path

        # Look for the date pattern in the filename
        parts = filename.split('_')
        if len(parts) >= 3:
            # The date should be in the last part before .txt or .tar.gz
            date_part = int(parts[2][:8])  # Extract YYYYMMDD

            # Add to dictionary
            if date_part in date_files_dict:
                date_files_dict[date_part].append(file)
            else:
                date_files_dict[date_part] = [file]

    return date_files_dict


@dg.sensor(
    job_name="unzip_archive_job",
)
def instrument_upload_sensor(context,
                             s3: S3ResourceNCSA
                             ):
    context.log.info("Instrument data upload")
    last_processed_date = int(context.cursor) if context.cursor else 0
    context.log.info(f"Last processed date: {last_processed_date}")
    files = group_files_by_date(
        list_files(EnvVar('DEST_BUCKET').get_value(), "raw", s3.get_client())
    )

    sensor_files = []
    sensor_date = None
    for data_date in sorted(files.keys()):
        if data_date <= last_processed_date:
            continue
        sensor_files = files[data_date]
        sensor_date = data_date
        break

    context.log.info(f"Found {len(sensor_files)} files on {sensor_date}")
    context.update_cursor(str(sensor_date))

    # After processing, there can be just the .txt file for that date
    if sensor_files and len(sensor_files) > 1:
        cloud_file = {}
        tar_gz_files = {}

        for file in sensor_files:
            filename = file.split('/')[-1]
            site_code = filename.split('_')[1]

            if filename.startswith("Cloud_"):
                cloud_file[site_code] = file
            elif "tar.gz" in file:
                if site_code in tar_gz_files:
                    tar_gz_files[site_code].append(file)
                else:
                    tar_gz_files[site_code] = [file]

        for site in tar_gz_files.keys():
            run_config = RunConfig({
                "unzip_chunked_archive": ChunkedArchiveConfig(
                    site=site,
                    observation_date=str(sensor_date),
                    cloud_file=cloud_file[site],
                    file_chunks=tar_gz_files[site]
                )
                }
            )
            print(f"{sensor_date}: {cloud_file[site]}, [{', '.join(tar_gz_files[site])}]")
            yield dg.RunRequest(run_key=f"sort-{sensor_date}-{site}", run_config=run_config)
