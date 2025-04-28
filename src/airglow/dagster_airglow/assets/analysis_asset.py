import os
import tempfile
from datetime import datetime, timedelta
from typing import Any

import dagster as dg

from dagster import EnvVar
from dagster_ncsa import S3ResourceNCSA
from dagster_mysql import MySQLResource


from airglow import FPIprocess, fpiinfo
from airglow.exceptions import NoLaserImagesError, NoSkyImagesError


def get_instrument_info(site, year, doy):
    """Get instrument information for the given site and date."""
    # Get date from year and doy
    nominal_dt = datetime(year, 1, 1) + timedelta(days=doy - 1)

    # Get the instrument name at this site
    instr_name = fpiinfo.get_instr_at(site, nominal_dt)[0]

    # Import the site information
    site_name = fpiinfo.get_site_of(instr_name, nominal_dt)

    # Create "minime05_uao_20130729" string
    datestr = nominal_dt.strftime('%Y%m%d')
    instrsitedate = instr_name + '_' + site_name + '_' + datestr

    return instr_name, site_name, datestr, instrsitedate


def download_fpi_data(context: dg.AssetExecutionContext,
                      s3: S3ResourceNCSA,
                      site, year, datestr, target_dir,
                      sky_line_tag: str,
                      fpi_data_path: str,
                      cloud_cover_path: str,
                      bucket_prefix_dir:str = ''):
    """Download FPI data files from cloud storage."""

    created_files = []

    # Get the instrument name for the site
    instr_name = get_instrument_info(site, year, int(datestr[-3:]))[0]

    # Download FPI data
    context.log.info(f"Downloading FPI data for {fpi_data_path}")

    files = s3.list_files(EnvVar("DEST_BUCKET").get_value(),
                          fpi_data_path, "")

    for f in files:
        if f"_{sky_line_tag}_" not in os.path.basename(f):
            continue

        relative_path = os.path.dirname(f[len(bucket_prefix_dir):])
        relative_file = f[len(bucket_prefix_dir):]

        local_path = os.path.join(target_dir, relative_path)
        os.makedirs(local_path, exist_ok=True)

        local_file = os.path.join(target_dir, relative_file)
        s3.get_client().download_file(
            Bucket=EnvVar("DEST_BUCKET").get_value(),
            Key=f,
            Filename=local_file
        )

        created_files.append(local_file)

    # Download cloud sensor data
    context.log.info(f"Downloading cloud sensor data for {site} on {datestr}")
    files = s3.list_files(
        EnvVar("DEST_BUCKET").get_value(),
        cloud_cover_path, "txt"
    )

    for f in files:
        relative_path = os.path.dirname(f[len(bucket_prefix_dir):])
        relative_file = f[len(bucket_prefix_dir):]

        local_path = os.path.join(target_dir, relative_path)
        os.makedirs(local_path, exist_ok=True)

        local_file = os.path.join(target_dir, relative_file)

        s3.get_client().download_file(
            Bucket=EnvVar("DEST_BUCKET").get_value(),
            Key=f,
            Filename=local_file
        )

        created_files.append(local_file)

    return created_files, instr_name


def upload_results(context: dg.AssetExecutionContext,
                   s3: S3ResourceNCSA,
                   local_dir: str, object_prefix: str):
    """Upload the results to cloud storage."""
    if not os.path.exists(local_dir):
        context.log.info(f"Local directory {local_dir} does not exist.")
        return

    # Upload the results to the cloud storage
    bucket_name = EnvVar("DEST_BUCKET").get_value()
    s3_client = s3.get_client()
    context.log.info(f"Uploading results to /{object_prefix}")

    # Ensure the prefix ends with a slash if it's not empty
    if object_prefix and not object_prefix.endswith('/'):
        object_prefix += '/'

    uploaded_files = []

    # Walk through the directory
    for root, dirs, files in os.walk(local_dir):
        for file in files:
            # Get the full local path
            local_path = os.path.join(root, file)

            # Calculate the relative path from the base directory
            relative_path = os.path.relpath(local_path, local_dir)

            # Create the S3 object key with the prefix
            s3_key = object_prefix + relative_path

            try:
                # Upload the file
                s3_client.upload_file(local_path, bucket_name, s3_key)
                uploaded_files.append(local_path)
                context.log.info(f"Uploaded {local_path} to s3://{bucket_name}/{s3_key}")
            except Exception as e:
                print(f"Error uploading {local_path}: {str(e)}")

    return uploaded_files


@dg.asset(
    name="analyze_data_x",
    ins={"unzip_chunked_archive": dg.AssetIn(dagster_type=dict[str, Any])},
)
def analyze_data_x(context: dg.AssetExecutionContext,
                   unzip_chunked_archive: dict[str, Any],
                   s3: S3ResourceNCSA,
                   mysql: MySQLResource) -> str:
    return analyze_data(context, unzip_chunked_archive, "X", s3, mysql)


@dg.asset(
    name="analyze_data_xr",
    ins={"unzip_chunked_archive": dg.AssetIn(dagster_type=dict[str, Any])},
)
def analyze_data_xr(context: dg.AssetExecutionContext,
                    unzip_chunked_archive: dict[str, Any],
                    s3: S3ResourceNCSA,
                    mysql: MySQLResource) -> str:
    return analyze_data(context, unzip_chunked_archive, "XR", s3, mysql)


@dg.asset(
    name="analyze_data_xg",
    ins={"unzip_chunked_archive": dg.AssetIn(dagster_type=dict[str, Any])},
)
def analyze_data_xg(context: dg.AssetExecutionContext,
                    unzip_chunked_archive: dict[str, Any],
                    s3: S3ResourceNCSA,
                    mysql: MySQLResource) -> str:
    return analyze_data(context, unzip_chunked_archive, "XG", s3, mysql)


def analyze_data(context: dg.AssetExecutionContext,
                 unzip_chunked_archive: dict[str, any],
                 sky_line_tag: str,
                 s3: S3ResourceNCSA,
                 mysql: MySQLResource) -> str:
    """
    Analyze the data and return the analysis result.
    :param context: The asset execution context.
    :param unzip_chunked_archive: Output from the unzip asset.
    :param s3: The S3 resource.
    :return: The analysis result.
    """
    # Perform some analysis on the data
    observation_date = unzip_chunked_archive["observation_date"]
    date_obj = datetime.strptime(observation_date, "%Y%m%d")
    doy = date_obj.timetuple().tm_yday
    year = int(unzip_chunked_archive["year"])

    # Get date string from year and day of year
    instr_name, site_name, datestr, instrsitedate = \
        get_instrument_info(
            unzip_chunked_archive["site"],
            year,
            doy
        )

    context.log.info("Instrument name: %s", instr_name)
    context.log.info("Site name: %s", site_name)
    context.log.info("Date string: %s", datestr)
    context.log.info("Instrument site date: %s", instrsitedate)

    # Create a temporary directory context manager
    with tempfile.TemporaryDirectory() as temp_dir:
        # Define all paths using the temporary directory
        fpi_dir = os.path.join(temp_dir, 'fpi/')
        bw_dir = os.path.join(temp_dir, 'cloudsensor/')
        x300_dir = os.path.join(temp_dir, 'templogs/x300/')
        results_stub = os.path.join(temp_dir, 'test/results/')
        madrigal_stub = os.path.join(temp_dir, 'test/madrigal/')
        share_stub = os.path.join(temp_dir, 'share/')
        temp_plots_stub = os.path.join(temp_dir, 'temporary_plots/')

        # Make sure all directories exist
        for directory in [fpi_dir, bw_dir, x300_dir, results_stub,
                          madrigal_stub, share_stub, temp_plots_stub]:
            os.makedirs(directory, exist_ok=True)

        download_fpi_data(
            context=context, 
            s3=s3, 
            site=unzip_chunked_archive["site"], 
            year=year,
            sky_line_tag=sky_line_tag,
            datestr=observation_date, 
            target_dir=temp_dir,
            fpi_data_path=unzip_chunked_archive['fpi_data_path'],
            cloud_cover_path=unzip_chunked_archive['cloud_cover_path']
        )

        try:
            FPIprocess.process_instr(instr_name, year, doy,
                                     mysql=mysql,
                                     sky_line_tag=sky_line_tag,
                                     fpi_dir=fpi_dir, bw_dir=bw_dir,
                                     send_to_madrigal=True,
                                     send_to_website=True,
                                     x300_dir=x300_dir, results_stub=results_stub,
                                     madrigal_stub=madrigal_stub, share_stub=share_stub,
                                     temp_plots_stub=temp_plots_stub)

            upload_results(context, s3, results_stub, EnvVar("RESULTS_PATH").get_value())
            upload_results(context, s3, temp_plots_stub,
                           EnvVar("SUMMARY_IMAGES_PATH").get_value())
            upload_results(context, s3, madrigal_stub,
                           EnvVar("MADRIGAL_PATH").get_value())
        except NoSkyImagesError:
            context.log.error(
                f"No SkyImages found for {sky_line_tag} on {observation_date}")
        except NoLaserImagesError:
            context.log.error(
                f"No LaserImages found for {sky_line_tag} on {observation_date}")

    return "ok"
