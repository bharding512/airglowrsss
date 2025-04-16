import os
import tempfile
from datetime import datetime, timedelta

import dagster as dg
from dagster import EnvVar
from dagster_ncsa import S3ResourceNCSA

from airglow import FPIprocess, fpiinfo


def get_instrument_info(site, year, doy):
    """Get instrument information for the given site and date."""
    # Get date from year and doy
    nominal_dt = datetime(year, 1, 1) + timedelta(days=doy - 1)

    # Get the instrument name at this site
    instr_name = fpiinfo.get_instr_at(site, nominal_dt)[0]

    # Import the site information
    site_name = fpiinfo.get_site_of(instr_name, nominal_dt)

    # Import the instrument information
    instrument = fpiinfo.get_instr_info(instr_name, nominal_dt)

    # Create "minime05_uao_20130729" string
    datestr = nominal_dt.strftime('%Y%m%d')
    instrsitedate = instr_name + '_' + site_name + '_' + datestr

    return instr_name, site_name, datestr, instrsitedate

class AnalysisConfig(dg.Config):
    site: str = "uao"
    year: int = 2025
    observation_date: str = "20250413"
    fpi_data_path: str = "fpi/minime05/uao/2025/20250413"
    cloud_cover_path: str = "cloudsensor/uao/2025"

def download_fpi_data(context: dg.AssetExecutionContext,
                      config: AnalysisConfig,
                      s3: S3ResourceNCSA,
                      site, year, datestr, target_dir,
                      bucket_prefix_dir=''):
    """Download FPI data files from cloud storage."""

    created_files = []

    # Get the instrument name for the site
    instr_name = get_instrument_info(site, year, int(datestr[-3:]))[0]

    # Download FPI data
    context.log.info(f"Downloading FPI data for {config.fpi_data_path}")

    files = s3.list_files(EnvVar("DEST_BUCKET").get_value(),
                          config.fpi_data_path, "hdf5")

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

    # Download cloud sensor data
    context.log.info(f"Downloading cloud sensor data for {site} on {datestr}")
    files = s3.list_files(
        EnvVar("DEST_BUCKET").get_value(),
        config.cloud_cover_path, "txt"
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

@dg.asset
def analyze_data(context: dg.AssetExecutionContext,
                 config: AnalysisConfig,
                 s3: S3ResourceNCSA) -> str:
    """
    Analyze the data and return the analysis result.
    :param context: The asset execution context.
    :param config: The data to analyze.
    :param s3: The S3 resource.
    :return: The analysis result.
    """
    # Perform some analysis on the data
    date_obj = datetime.strptime(config.observation_date, "%Y%m%d")
    doy = date_obj.timetuple().tm_yday

    # Get date string from year and day of year
    instr_name, site_name, datestr, instrsitedate = get_instrument_info(config.site,
                                                                        config.year, doy)

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
            context, config, s3, config.site, config.year,
            config.observation_date, temp_dir, bucket_prefix_dir=""
        )

        FPIprocess.process_instr(instr_name, config.year, doy,
                                 sky_line_tag="XR",
                                 fpi_dir=fpi_dir, bw_dir=bw_dir,
                                 send_to_madrigal=True,
                                 x300_dir=x300_dir, results_stub=results_stub,
                                 madrigal_stub=madrigal_stub, share_stub=share_stub,
                                 temp_plots_stub=temp_plots_stub)

        upload_results(context, s3, results_stub, EnvVar("RESULTS_PATH").get_value())
        upload_results(context, s3, temp_plots_stub, EnvVar("SUMMARY_IMAGES_PATH").get_value())
        upload_results(context, s3, madrigal_stub, EnvVar("MADRIGAL_PATH").get_value())
    return "ok"