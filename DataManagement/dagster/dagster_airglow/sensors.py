from pathlib import Path

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
    """
    Using a list of files from the S3 bucket, this function groups the files by date.
    """
    date_files_dict = {}

    for file in file_list:
        # Extract the date part from the filename
        # The date appears in format YYYYMMDD
        filename = Path(file).name

        # Look for the date pattern in the filename
        parts = filename.split('_')
        if len(parts) >= 3:
            # The date should be in the last part before .txt or .tar.gz
            date_part = int(parts[2][:8])  # Extract YYYYMMDD

            # Add to dictionary or append to existing entry if there
            if date_part in date_files_dict:
                date_files_dict[date_part].append(file)
            else:
                date_files_dict[date_part] = [file]

    return date_files_dict


def cloud_cover_files_for_site(site: str, files: list[str]) -> list[str]:
    """
    Filters the list of files to get only the cloud cover files for the specified site.
    """
    return [file for file in files if f"Cloud_{site}" in file and file.endswith(".txt")]


@dg.sensor(
    job_name="unzip_archive_job",
    minimum_interval_seconds=1 * 60 * 60,  # 1 hour
)
def instrument_upload_sensor(context,
                             s3: S3ResourceNCSA
                             ):
    objects = list_files(EnvVar('DEST_BUCKET').get_value(), "raw", s3.get_client())
    files = group_files_by_date(objects)

    for data_date in sorted(files.keys()):
        sensor_files = files[data_date]
        sensor_date = data_date

        context.log.info(f"Found {len(sensor_files)} files on {sensor_date}")

        # After processing, there can be just the .txt file for that date
        if sensor_files and len(sensor_files) > 1:
            tar_gz_files = {}
            complete_sites = set()

            for file in sensor_files:
                filename = file.split('/')[-1]
                site_code = filename.split('_')[1]

                if filename.startswith("fpi05") and filename.endswith(".txt"):
                    # This is a log file, we can ignore it, but it signifies that all
                    # the files are uploaded for this date/site
                    complete_sites.add(site_code)
                    continue

                if "tar.gz" in file:
                    if site_code in tar_gz_files:
                        tar_gz_files[site_code].append(file)
                    else:
                        tar_gz_files[site_code] = [file]

            for site in tar_gz_files.keys():
                if site in complete_sites:
                    run_config = RunConfig({
                        "unzip_chunked_archive": ChunkedArchiveConfig(
                            site=site,
                            observation_date=str(sensor_date),
                            cloud_files=cloud_cover_files_for_site(site, objects),
                            file_chunks=tar_gz_files[site]
                        )
                        }
                    )
                    yield dg.RunRequest(
                        run_key=f"sort-{sensor_date}-{site}",
                        run_config=run_config
                    )
                else:
                    context.log.info(f"Incomplete upload for {site} on {sensor_date} - will pick them up next time")  # NOQA E501
