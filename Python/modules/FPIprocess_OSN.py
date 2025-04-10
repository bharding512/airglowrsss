#!/usr/bin/env python3
# FPIprocess_OSN.py
# Script to download FPI data from OSN, process it, and upload results

import os
import sys
import argparse
import datetime as datetime
import logging
from pathlib import Path
import fpiinfo
import FPIprocess
from cloud_storage import CloudStorage, Configuration


def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


def download_fpi_data(cloud_storage, config, site, year, datestr, target_dir, bucket_prefix_dir=''):
    """Download FPI data files from cloud storage."""
    logger = logging.getLogger(__name__)
    created_files = []
    
    # Get the instrument name for the site
    instr_name = get_instrument_info(site, year, int(datestr[-3:]))[0]
    
    # Download FPI data
    logger.info(f"Downloading FPI data for {instr_name} at {site} on {datestr}")
    fpi_prefix = f'{config.aws_fpi_prefix}{instr_name}/{site}/{year}/{datestr}/'
    files = cloud_storage.list_objects(fpi_prefix)
    
    for f in files:
        relative_path = os.path.dirname(f[len(bucket_prefix_dir):])
        relative_file = f[len(bucket_prefix_dir):]

        local_path = os.path.join(target_dir, relative_path)
        os.makedirs(local_path, exist_ok=True)

        local_file = os.path.join(target_dir, relative_file)

        cloud_storage.download_file(f, local_file)
        created_files.append(local_file)

    # Download cloud sensor data
    logger.info(f"Downloading cloud sensor data for {site} on {datestr}")
    cloud_prefix = f'{config.aws_cloudsensor_prefix}{site}/{year}/Cloud_{site}'
    files = cloud_storage.list_objects(cloud_prefix)
    
    for f in files:
        relative_path = os.path.dirname(f[len(bucket_prefix_dir):])
        relative_file = f[len(bucket_prefix_dir):]

        local_path = os.path.join(target_dir, relative_path)
        os.makedirs(local_path, exist_ok=True)

        local_file = os.path.join(target_dir, relative_file)

        cloud_storage.download_file(f, local_file)
        created_files.append(local_file)
    
    return created_files, instr_name


def delete_empty_directories(root_dir):
    """Delete empty directories recursively."""
    logger = logging.getLogger(__name__)
    removed_dirs = False
    for dirpath, dirnames, filenames in os.walk(root_dir, topdown=False):
        if not dirnames and not filenames:
            os.rmdir(dirpath)
            logger.info(f"Deleted empty directory: {dirpath}")
            removed_dirs = True

    return removed_dirs


def delete_files_and_empty_dirs(file_list, target_dir):
    """Delete all files in the provided list and remove empty directories."""
    logger = logging.getLogger(__name__)
    
    # Track directories to check later
    dirs_to_check = set()
    
    # First delete all files
    for file_path in file_list:
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
                logger.info(f"Successfully deleted: {file_path}")
                
                # Add the directory to our check list
                dir_path = os.path.dirname(file_path)
                dirs_to_check.add(dir_path)
            else:
                logger.warning(f"File not found: {file_path}")
        except Exception as e:
            logger.error(f"Error deleting {file_path}: {e}")

    # Recursively delete empty directories
    result = delete_empty_directories(target_dir)
    while result:
        result = delete_empty_directories(target_dir)


def upload_directory_to_cloud(local_dir, cloud_prefix, cloud_storage):
    """Upload a directory and its contents to cloud storage. Delete local files after successful upload."""
    logger = logging.getLogger(__name__)
    local_dir = Path(local_dir)
    
    if not local_dir.exists():
        logger.warning(f"Directory does not exist: {local_dir}")
        return
        
    for root, _, files in os.walk(local_dir):
        for file in files:
            local_path = Path(root) / file
            relative_path = local_path.relative_to(local_dir)
            cloud_key = f"{cloud_prefix}{relative_path}"
            
            try:
                # Upload the file
                logger.info(f"Uploading {str(local_path)}")
                cloud_storage.upload_file(str(local_path), cloud_key)
                
                # If upload is successful, delete the local file
                local_path.unlink()
                
            except Exception as e:
                # Log error and continue with other files if one fails
                logger.error(f"Error uploading {local_path}: {str(e)}")


def clean_up_post_FPIprocess(created_files, target_dir, config, cloud_storage):
    """Clean up temporary files and upload results to cloud storage."""
    logger = logging.getLogger(__name__)
    
    # Delete all FPI data files
    delete_files_and_empty_dirs(created_files, target_dir)

    # Upload processed results to cloud storage
    results_dir = Path(config.fpi_process_kwargs['results_stub'])
    if results_dir.exists():
        cloud_results_prefix = f"{config.aws_results_prefix}"
        upload_directory_to_cloud(results_dir, cloud_results_prefix, cloud_storage)
        logger.info(f"Uploaded results to cloud storage: {cloud_results_prefix}")

    # Upload temporary plots to cloud storage
    temp_plots_dir = Path(config.fpi_process_kwargs['temp_plots_stub'])
    if temp_plots_dir.exists():
        cloud_results_prefix = f"{config.aws_summaryimages_prefix}"
        upload_directory_to_cloud(temp_plots_dir, cloud_results_prefix, cloud_storage)
        logger.info(f"Uploaded plots to cloud storage: {cloud_results_prefix}")

    # Upload madrigal files to cloud storage
    madrigal_dir = Path(config.fpi_process_kwargs['madrigal_stub'])
    if madrigal_dir.exists():
        cloud_results_prefix = f"{config.aws_madrigal_prefix}"
        upload_directory_to_cloud(madrigal_dir, cloud_results_prefix, cloud_storage)
        logger.info(f"Uploaded madrigal files to cloud storage: {cloud_results_prefix}")


def get_instrument_info(site, year, doy):
    """Get instrument information for the given site and date."""
    # Get date from year and doy
    nominal_dt = datetime.datetime(year, 1, 1) + datetime.timedelta(days=doy-1)
    
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


def process_fpi_data(site, year, doy, env_file='.env'):
    """Main function to process FPI data."""
    logger = setup_logging()
    
    # Load configuration
    config = Configuration(env_file)
    cloud_storage = CloudStorage(config)
    
    # Get date string from year and day of year
    instr_name, site_name, datestr, instrsitedate = get_instrument_info(site, year, doy)
    
    logger.info(f"Processing {instrsitedate}")
    
    # Set up target directory
    target_dir = config.temp_dir
    bucket_prefix_dir = ''
    
    # Download data
    created_files, instr_name = download_fpi_data(
        cloud_storage, config, site, year, datestr, target_dir, bucket_prefix_dir
    )
    
    # Process the FPI data
    logger.info(f"Processing FPI data for {instr_name} on day {doy} of {year}")
    FPIprocess.process_instr(instr_name, year, doy, **config.fpi_process_kwargs)
    
    # Clean up and upload results
    logger.info("Cleaning up and uploading results")
    clean_up_post_FPIprocess(created_files, target_dir, config, cloud_storage)
    
    logger.info(f"FPI processing complete for {instrsitedate}")


def main():
    """Parse command line arguments and process FPI data."""
    parser = argparse.ArgumentParser(description='Process FPI data from OSN cloud storage.')
    parser.add_argument('-s', '--site', required=True, help='Site identifier (e.g., uao)')
    parser.add_argument('-y', '--year', required=True, type=int, help='Year (e.g., 2025)')
    parser.add_argument('-d', '--doy', required=True, type=int, help='Day of year (e.g., 99)')
    parser.add_argument('-e', '--env', default='.env', help='Path to .env file')
    
    args = parser.parse_args()
    
    process_fpi_data(args.site, args.year, args.doy, args.env)


if __name__ == "__main__":
    main()
