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
    instrument_path: str = "fpi/minime05/uao/2025/20250413"
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

    files = s3.list_files(
        EnvVar("DEST_BUCKET").get_value(),
        config.instrument_path, extension="hdf5")
    context.log.info("Found %s files", list(files))

    # Create a temporary directory context manager
    with tempfile.TemporaryDirectory() as temp_dir:
        # Define all paths using the temporary directory
        fpi_dir = os.path.join(temp_dir, 'mango/fpi/')
        bw_dir = os.path.join(temp_dir, 'mango/templogs/cloudsensor/')
        x300_dir = os.path.join(temp_dir, 'mango/templogs/x300/')
        results_stub = os.path.join(temp_dir, 'mango/results/')
        madrigal_stub = os.path.join(temp_dir, 'mango/madrigal/')
        share_stub = os.path.join(temp_dir, 'mango/share/')
        temp_plots_stub = os.path.join(temp_dir, 'mango/temporary_plots/')

        # Make sure all directories exist
        for directory in [fpi_dir, bw_dir, x300_dir, results_stub,
                          madrigal_stub, share_stub, temp_plots_stub]:
            os.makedirs(directory, exist_ok=True)

    FPIprocess.process_instr(instr_name, config.year, doy,
                             fpi_dir=fpi_dir, bw_dir=bw_dir,
                             x300_dir=x300_dir, results_stub=results_stub,
                             madrigal_stub=madrigal_stub, share_stub=share_stub,
                             temp_plots_stub=temp_plots_stub)

    return "ok"