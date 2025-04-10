#!/usr/bin/env python3
"""
MANGO Data Pipeline Sorter (Cloud Version)

This module handles the sorting and processing of instrument data files,
primarily focused on FPI and Cloud sensor data with extensibility for other instruments.
Now works with AWS cloud storage.
"""

import os
import sys
import shutil
import tarfile
import tempfile
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple
from dataclasses import dataclass
import logging
from abc import ABC, abstractmethod
import numpy as np
import pytz
import boto3
from botocore.exceptions import ClientError
from botocore.config import Config
import io
import re

# Third-party imports
from optparse import OptionParser

# Project imports
from Zipper import activeinstruments
import fpiinfo
import FPI
import FPIprocess
from GitHubIssueHandler import SiteIssueManager, IssueType
from github_config import github_config
from exceptions import *

from dotenv import load_dotenv

def setup_logging(site: str):
    """Configure logging for the Sorter."""
    # Create logs directory if it doesn't exist
    log_dir = '/home/jmakela/logs'
    os.makedirs(log_dir, exist_ok=True)
    
    # Add filter to only allow __main__ logger
    class MainLoggingFilter(logging.Filter):
        def filter(self, record):
            allowed_loggers = {'__main__', 'py.warnings', 'root'}
            return record.name in allowed_loggers
    
    # Configure file handler with rotation
    log_file = f"{log_dir}/sorter_{site}.log"
    from logging.handlers import RotatingFileHandler
    file_handler = RotatingFileHandler(
        log_file,
        maxBytes=10*1024*1024,  # 10MB
        backupCount=5
    )
    file_handler.addFilter(MainLoggingFilter())
    
    console_handler = logging.StreamHandler()
    console_handler.addFilter(MainLoggingFilter())
    
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    logger.info(f"Starting Sorter for site: {site}")
    return logger

@dataclass
class FileMetadata:
    """Represents metadata extracted from an info file."""
    instrument: str
    site: str
    date: datetime
    instrument_number: str
    original_filename: str
    expected_parts: int
    total_size: int
    emails: List[str]

class ProcessingError(Exception):
    """Base class for processing errors with context."""
    def __init__(self, message: str, context: dict):
        super().__init__(message)
        self.context = context
        self.warning_log = []

class CloudStorage:
    """Interface for AWS cloud storage operations."""
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Initialize S3v4 client
        s3v4_config = Config(
            signature_version='s3v4',
            s3={'addressing_style': 'path'}
        )

        s3v4_client_args = {
            'service_name': 's3',
            'aws_access_key_id': config.aws_access_key,
            'aws_secret_access_key': config.aws_secret_key,
            'config': s3v4_config,
            'endpoint_url': config.aws_endpoint_url
        }

        self.client_s3v4 = boto3.client(**s3v4_client_args)

        # Initialize S3 client
        s3_config = Config(
            signature_version='s3',
            s3={'addressing_style': 'path'}
        )

        s3_client_args = {
            'service_name': 's3',
            'aws_access_key_id': config.aws_access_key,
            'aws_secret_access_key': config.aws_secret_key,
            'config': s3_config,
            'endpoint_url': config.aws_endpoint_url
        }

        self.client_s3 = boto3.client(**s3_client_args)

        self.bucket = config.dest_bucket
        
    def list_objects(self, prefix: str) -> List[str]:
        """List objects in the bucket with the given prefix."""
        try:
            response = self.client_s3v4.list_objects_v2(
                Bucket=self.bucket,
                Prefix=prefix
            )
            
            if 'Contents' in response:
                return [obj['Key'] for obj in response['Contents']]
            return []
            
        except ClientError as e:
            self.logger.error(f"Error listing objects with prefix {prefix}: {str(e)}")
            return []
    
    def download_file(self, key: str, local_path: str) -> bool:
        """Download file from cloud storage to local path."""
        try:
            self.logger.info(f"Downloading {key} to {local_path}")
            os.makedirs(os.path.dirname(local_path), exist_ok=True)
            self.client_s3.download_file(self.bucket, key, local_path)
            return True
        except ClientError as e:
            self.logger.error(f"Error downloading {key}: {str(e)}")
            return False
    
    def upload_file(self, local_path: str, key: str) -> bool:
        """Upload file from local path to cloud storage."""
        try:
            self.logger.info(f"Uploading {local_path} to {self.bucket},{key}")
            self.client_s3.upload_file(local_path, self.bucket, key)
            return True
        except ClientError as e:
            self.logger.error(f"Error uploading {local_path}: {str(e)}")
            return False
    
    def download_fileobj(self, key: str) -> Optional[io.BytesIO]:
        """Download file object from cloud storage."""
        try:
            fileobj = io.BytesIO()
            self.client_s3.download_fileobj(self.bucket, key, fileobj)
            fileobj.seek(0)
            return fileobj
        except ClientError as e:
            self.logger.error(f"Error downloading fileobj {key}: {str(e)}")
            return None
    
    def upload_fileobj(self, fileobj: io.BytesIO, key: str) -> bool:
        """Upload file object to cloud storage."""
        try:
            fileobj.seek(0)
            self.client_s3.upload_fileobj(fileobj, self.bucket, key)
            return True
        except ClientError as e:
            self.logger.error(f"Error uploading fileobj to {key}: {str(e)}")
            return False
    
    def read_text_file(self, key: str) -> Optional[str]:
        """Read text file from cloud storage."""
        fileobj = self.download_fileobj(key)
        if fileobj:
            return fileobj.read().decode('utf-8')
        return None
    
    def delete_object(self, key: str) -> bool:
        """Delete object from cloud storage."""
        try:
            self.client_s3.delete_object(Bucket=self.bucket, Key=key)
            return True
        except ClientError as e:
            self.logger.error(f"Error deleting {key}: {str(e)}")
            return False

class Configuration:
    """Configuration management."""
    def __init__(self):
        # AWS Cloud Configuration
        self.aws_access_key = os.getenv('AWS_ACCESS_KEY')
        self.aws_secret_key = os.getenv('AWS_SECRET_KEY')
        self.aws_endpoint_url = os.getenv('AWS_ENDPOINT_URL', 'https://ncsa.osn.xsede.org')
        self.dest_bucket = os.getenv('DEST_BUCKET', 'airglow')
        
        # Cloud directory structure
        self.aws_rx_prefix = 'raw/'       # was airglow/raw
        self.aws_data_prefix = '' # was airglow/
        self.aws_results_prefix = 'results/' # was airglow/fpi/results
        self.aws_summaryimages_prefix = 'SummaryImages/'
        self.aws_madrigal_prefix = 'madrigal/'
        
        # Local temporary directory
        self.temp_dir = os.getenv('TEMP_DIR', '/home/jmakela/tmp/mango')
        os.makedirs(self.temp_dir, exist_ok=True)
        
        # Processing configurations
        self.fpi_process_kwargs = {
            'reference': 'laser',
            'send_to_website': True,
            'enable_share': False,
            'send_to_madrigal': True,
            'sky_line_tag': 'X',
            'fpi_dir': f"{self.temp_dir}/fpi/",
            'bw_dir': f"{self.temp_dir}/cloudsensor/", # was /templogs/cloudsensor
            'x300_dir': f"{self.temp_dir}/templogs/x300/",
            'results_stub': f"{self.temp_dir}/results/",
            'temporary_plots_stub': f"{self.temp_dir}/temporary_plots/",
            'madrigal_stub': f"{self.temp_dir}/madrigal/",
        }

class InstrumentHandler(ABC):
    """Base class for instrument-specific handlers."""
    def __init__(self, config: Configuration, cloud_storage: CloudStorage):
        self.config = config
        self.cloud_storage = cloud_storage
        self.logger = logging.getLogger(__name__)
        
    @abstractmethod
    def sort_files(self, metadata: FileMetadata) -> List[str]:
        """Sort files for this instrument type."""
        pass
        
    @abstractmethod
    def process_data(self, metadata: FileMetadata, sorted_files: List[str]):
        """Process the sorted data."""
        pass

    def get_cloud_target_prefix(self, metadata: FileMetadata) -> str:
        """Get the target prefix for the instrument data in cloud storage."""
        return self.config.aws_data_prefix

    def get_local_target_directory(self, metadata: FileMetadata) -> Path:
        """Get the local target directory for the instrument data."""
        return Path(self.config.temp_dir)

class FPIHandler(InstrumentHandler):
    """Handler for FPI instruments."""
    def __init__(self, config: Config, cloud_storage: CloudStorage):
        super().__init__(config, cloud_storage)
        self.issue_manager = SiteIssueManager(
            token=github_config.github_token,
            repo_name=github_config.github_repo
        )

    def get_cloud_target_prefix(self, metadata: FileMetadata) -> str:
        """Get the target prefix for FPI data in cloud storage."""
        return f"{self.config.aws_data_prefix}fpi/minime{metadata.instrument_number}/{metadata.site}/{metadata.date.year}/"

    def get_local_target_directory(self, metadata: FileMetadata) -> Path:
        """Get the local target directory for FPI data."""
        return Path(self.config.temp_dir) / 'fpi' / f"minime{metadata.instrument_number}" / metadata.site / str(metadata.date.year)

    def sort_files(self, metadata: FileMetadata) -> List[str]:
        target_dir = self.get_local_target_directory(metadata)
        self.logger.info(f"Sorting files to local directory: {target_dir}")
        
        target_dir.mkdir(parents=True, exist_ok=True)
        result = self._concat_and_extract(metadata, target_dir)
        
        if result:
            self.logger.info(f"Successfully sorted {len(result)} files")
            self._set_permissions(target_dir, result)
            
            # Upload sorted files to cloud storage
            cloud_prefix = self.get_cloud_target_prefix(metadata)
            for file in result:
                local_file = target_dir / file
                cloud_key = f"{cloud_prefix}{file}"
                self.cloud_storage.upload_file(str(local_file), cloud_key)
                
            self.logger.info(f"Uploaded {len(result)} files to cloud storage")
        else:
            self.logger.warning("No files were sorted")
        
        return result

    def process_data(self, metadata: FileMetadata, sorted_files: List[str]):
        if not sorted_files:
            self.logger.warning("No files to process")
            return
            
        self.current_metadata = metadata  # Store for path resolution
        try:
            ldn = self._get_local_datetime(sorted_files[0])
            year, doy = self._calculate_processing_date(ldn, metadata.site)
            self.logger.info(f"Processing data for year={year}, doy={doy}")

            for tag, kwarg_tag in [('XG_', 'XG'), ('XR_', 'XR'), ('X_', 'X')]:
                if self._has_tag(sorted_files, tag):
                    self.logger.info(f"Processing tag: {tag}")
                    try:
                        self._process_with_tag(metadata, year, doy, kwarg_tag)
                        self.logger.info(f"Successfully processed {tag}")
                    except Exception as e:
                        self.logger.error(f"Error processing {tag}: {str(e)}", exc_info=True)
                        self._handle_processing_error(e, metadata, year, doy, kwarg_tag)
                        
            # Upload processed results to cloud storage
            results_dir = Path(self.config.fpi_process_kwargs['results_stub'])
            if results_dir.exists():
                cloud_results_prefix = f"{self.config.aws_results_prefix}/"
                self._upload_directory_to_cloud(results_dir, cloud_results_prefix)
                self.logger.info(f"Uploaded results to cloud storage: {cloud_results_prefix}")

            # Upload temporary plots to cloud storage
            temp_plots_dir = Path(self.config.fpi_process_kwargs['temporary_plots_stub'])
            if temp_plots_dir.exists():
                cloud_results_prefix = f"{self.config.aws_summaryimages_prefix}/"
                self._upload_directory_to_cloud(temp_plots_dir, cloud_results_prefix)
                self.logger.info(f"Uploaded plots to cloud storage: {cloud_results_prefix}")

            # Upload madrigal files to cloud storage
            madrigal_dir = PATH(self.config.fpi_process_kwargs['madrigal_stub'])
            if madrigal_dir.exists():
                cloud_results_prefix = f"{self.config.aws_madrigal_prefix}/"
                self._upload_directory_to_cloud(madrigal_dir, cloud_results_prefix)
                self.logger.info(f"Uploaded madrigal files to cloud storage: {cloud_results_prefix}")
        finally:
            self.current_metadata = None
            
    def upload_directory_to_cloud(self, local_dir: Path, cloud_prefix: str):
        """Upload a directory and its contents to cloud storage. Delete local files after successful upload."""
        for root, _, files in os.walk(local_dir):
            for file in files:
                local_path = Path(root) / file
                relative_path = local_path.relative_to(local_dir)
                cloud_key = f"{cloud_prefix}{relative_path}"
                
                try:
                    # Upload the file
                    self.cloud_storage.upload_file(str(local_path), cloud_key)
                    
                    # If upload is successful, delete the local file
                    local_path.unlink()
                    
                except Exception as e:
                    # Log error and continue with other files if one fails
                    self.logger.error(f"Error uploading {local_path}: {str(e)}")

    def _process_with_tag(self, metadata: FileMetadata, year: int, doy: int, tag: str):
        """Process data with a specific emission line tag."""
        process_kwargs = self.config.fpi_process_kwargs.copy()
        process_kwargs['sky_line_tag'] = tag
        
        warning = FPIprocess.process_instr(
            f"minime{metadata.instrument_number}", 
            year, 
            doy, 
            **process_kwargs
        )
        
        if warning:
            self._handle_warning(warning, metadata, year, doy, tag)

    def _handle_warning(self, warning: str, metadata: FileMetadata, year: int, doy: int, tag: str):
        """Handle processing warnings by creating GitHub issues."""
        self.issue_manager.handle_processing_issue(
            site_id=metadata.site,
            message=warning,
            category=IssueType.WARNING,
            additional_context={
                "Instrument": f"minime{metadata.instrument_number}",
                "Year": year,
                "Day of Year": doy,
                "Tag": tag
            }
        )

    def _concat_and_extract(self, metadata: FileMetadata, target_dir: Path) -> List[str]:
        """Concatenate split files from cloud storage and extract contents."""
        file_pattern = f"{metadata.original_filename}"
        rx_prefix = self.config.aws_rx_prefix
        
        print(f"concat_and_extract: {file_pattern}")

        # List all parts in cloud storage
        parts = self.cloud_storage.list_objects(f"{rx_prefix}{file_pattern}")
        parts = [p for p in parts if p.startswith(f"{rx_prefix}{metadata.original_filename}")]
        
        if len(parts) != metadata.expected_parts:
            self.logger.warning(f"Expected {metadata.expected_parts} parts, found {len(parts)}")
            return []

        temp_file = target_dir / f"temp_{datetime.now().strftime('%Y%m%d%H%M%S%f')}.tar.gz"
        try:
            # Concatenate parts
            with open(temp_file, 'wb') as outfile:
                for part in sorted(parts):
                    file_obj = self.cloud_storage.download_fileobj(part)
                    if file_obj:
                        outfile.write(file_obj.read())
                    else:
                        self.logger.error(f"Failed to download part: {part}")
                        return []

            # Verify size
            if temp_file.stat().st_size != metadata.total_size:
                self.logger.error(f"Size mismatch after concatenation: expected {metadata.total_size}, got {temp_file.stat().st_size}")
                return []

            # Extract
            with tarfile.open(temp_file, 'r:gz') as tar:
                tar.extractall(path=target_dir)
                extracted_files = [member.name for member in tar.getmembers()]

            # Clean up parts after successful extraction
            for part in parts:
                self.cloud_storage.delete_object(part)
            self.logger.info(f"Cleaned up {len(parts)} partial files from cloud storage")

            return extracted_files

        finally:
            if temp_file.exists():
                temp_file.unlink()

    def _get_local_datetime(self, first_file: str) -> datetime:
        """Get local datetime from first file."""
        full_path = self.get_local_target_directory(self.current_metadata) / first_file
        img_info = FPI.ReadIMG(str(full_path)).info
        return img_info['LocalTime']

    def _calculate_processing_date(self, ldn: datetime, site: str) -> tuple[int, int]:
        """Calculate processing year and day of year."""
        site_info = fpiinfo.get_site_info(site, ldn)
        utdn = ldn.replace(tzinfo=pytz.timezone(site_info['Timezone'])).astimezone(pytz.utc).replace(tzinfo=None)
        site_lon = np.mod(site_info['Location'][1]+180, 360)-180
        sltdn = utdn + timedelta(hours=24*site_lon/360.)
        dn0 = sltdn - timedelta(hours=12)
        return dn0.year, dn0.timetuple().tm_yday

    def _has_tag(self, files: List[str], tag: str) -> bool:
        """Check if any files have the given tag."""
        return any(tag in Path(f).name for f in files)

    def _set_permissions(self, target_dir: Path, files: List[str]):
        """Set appropriate permissions on processed files."""
        for file in files:
            file_path = target_dir / file
            os.chmod(file_path, 0o644)
            # Skip chown in cloud environment as it might not be applicable
            # or use a different approach if needed

    def _handle_processing_error(self, error: Exception, metadata: FileMetadata, 
                               year: int, doy: int, tag: str):
        """Handle processing errors."""
        self.issue_manager.handle_processing_issue(
            site_id=metadata.site,
            message=str(error),
            category=IssueType.ERROR,
            additional_context={
                "Instrument": f"minime{metadata.instrument_number}",
                "Year": year,
                "Day of Year": doy,
                "Tag": tag,
                "Error Type": type(error).__name__
            }
        )

class CloudSensorHandler(InstrumentHandler):
    """Handler for cloud sensor data."""
    def get_cloud_target_prefix(self, metadata: FileMetadata) -> str:
        """Get the target prefix for cloud sensor data in cloud storage."""
        return f"{self.config.aws_data_prefix}cloudsensor/{metadata.site}/" # was templogs/cloudsensor

    def get_local_target_directory(self, metadata: FileMetadata) -> Path:
        """Get the local target directory for cloud sensor data."""
        return Path(self.config.temp_dir) / 'cloudsensor' / metadata.site # was templogs/cloudsensor

    def sort_files(self, metadata: FileMetadata) -> List[str]:
        target_dir = self.get_local_target_directory(metadata)
        target_dir.mkdir(parents=True, exist_ok=True)
        
        # Download from cloud storage
        source_key = f"{self.config.aws_rx_prefix}{metadata.original_filename}"
        target_file = target_dir / metadata.original_filename
        
        if self.cloud_storage.download_file(source_key, str(target_file)):
            # Upload to the archive location
            target_key = f"{self.get_cloud_target_prefix(metadata)}{metadata.original_filename}"
            self.cloud_storage.upload_file(str(target_file), target_key)
            self.logger.info(f"Uploaded cloud sensor file to {target_key}")
            
            return [str(target_file)]
        
        return []

    def process_data(self, metadata: FileMetadata, sorted_files: List[str]):
        """Process cloud sensor data - implementation as needed."""
        # For now, just moving the file is sufficient which happens in sort_files
        pass

    def process_cloud_file(self, cloud_file_key: str):
        """Process a single cloud sensor file from cloud storage."""
        # Extract site from filename (Cloud_site_date.txt)
        filename = os.path.basename(cloud_file_key)
        match = re.match(r'Cloud_([^_]+)_.*\.txt', filename)
        
        if not match:
            self.logger.warning(f"Cannot extract site from cloud sensor filename: {filename}")
            return
            
        site = match.group(1)
        
        # Create target directories
        local_target_dir = Path(self.config.temp_dir) / 'cloudsensor' / site # was templogs/cloudsensor
        local_target_dir.mkdir(parents=True, exist_ok=True)
        
        local_file = local_target_dir / filename
        
        # Download from cloud storage
        if self.cloud_storage.download_file(cloud_file_key, str(local_file)):
            # Upload to the proper location in cloud storage
            cloud_target_key = f"{self.config.aws_data_prefix}cloudsensor/{site}/{filename}" # was templogs/cloudsensor
            if self.cloud_storage.upload_file(str(local_file), cloud_target_key):
                # Delete original file from rx directory
                self.cloud_storage.delete_object(cloud_file_key)
                self.logger.info(f"Processed cloud sensor file: {filename} for site {site}")
            else:
                self.logger.error(f"Failed to upload cloud sensor file to target location: {cloud_target_key}")
        else:
            self.logger.error(f"Failed to download cloud sensor file: {cloud_file_key}")

class DataProcessor:
    """Main data processing orchestrator."""
    def __init__(self, config: Configuration):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.cloud_storage = CloudStorage(config)
        self.handlers = {
            'fpi': FPIHandler(config, self.cloud_storage),
            'cloud': CloudSensorHandler(config, self.cloud_storage),
        }

    def process_site(self, site: str):
        """Process all files for a given site."""
        try:
            self._acquire_lock(site)
            self.logger.info(f"Processing site: {site}")
            
            # Process cloud sensor files first
            cloud_files = self.get_cloud_files(site)
            if cloud_files:
                self.logger.info(f"Found {len(cloud_files)} cloud sensor files")
                cloud_handler = self.handlers['cloud']
                for cloud_file in cloud_files:
                    try:
                        cloud_handler.process_cloud_file(cloud_file)
                    except Exception as e:
                        self.logger.error(f"Error processing cloud file {cloud_file}: {str(e)}", exc_info=True)
            
            # Process other instrument files
            info_files = self._get_info_files(site)
            self.logger.info(f"Found {len(info_files)} info files to process")
            
            for info_file in info_files:
                try:
                    self.logger.info(f"Processing file: {info_file}")
                    metadata = self._parse_info_file(info_file)
#                    metadata = FileMetadata(
#                        instrument='fpi',
#                        site=site,
#                        date = datetime.strptime('2025-03-25 02:46:09', '%Y-%m-%d %H:%M:%S'),
#                        instrument_number='05',
#                        original_filename=None,
#                        expected_parts=None,
#                        total_size=None,
#                        emails=None
#                    )
                    handler = self.handlers.get(metadata.instrument)
                    if handler:
                        self.logger.info(f"Using handler for instrument: {metadata.instrument}")
                        sorted_files = handler.sort_files(metadata)
#                        import glob as glob
#                        sorted_files = glob.glob('/tmp/mango/fpi/minime05/uao/2025/20250325/*hdf5')
#                        print(sorted_files)
                        if sorted_files:
                            self.logger.info(f"Sorted {len(sorted_files)} files")
                            # Delete info file from cloud storage
                            self.cloud_storage.delete_object(info_file)
                            handler.process_data(metadata, sorted_files)
                            self.logger.info("Processing complete")
                            
                            # Clean up local temporary files
                            self._cleanup_temp_files(handler, metadata)
                    else:
                        self.logger.warning(f"No handler for instrument: {metadata.instrument}")
                except Exception as e:
                    self.logger.error(f"Error processing {info_file}: {str(e)}", exc_info=True)
                    self._handle_error(e, info_file)
        finally:
            self._release_lock(site)
            self.logger.info(f"Completed processing for site: {site}")

    def _cleanup_temp_files(self, handler, metadata):
        """Clean up local temporary files after processing."""
        temp_dir = handler.get_local_target_directory(metadata)
        if temp_dir.exists() and temp_dir.is_dir():
            self.logger.info(f"Cleaning up temporary directory: {temp_dir}")
            try:
                shutil.rmtree(str(temp_dir))
            except Exception as e:
                self.logger.warning(f"Error cleaning up directory {temp_dir}: {str(e)}")

    def get_cloud_files(self, site: str) -> List[str]:
        """Get cloud sensor files for a site from cloud storage."""
        rx_prefix = self.config.aws_rx_prefix
        if site == 'all':
            pattern = f"{rx_prefix}Cloud_"
        else:
            pattern = f"{rx_prefix}Cloud_{site}_"
            
        all_files = self.cloud_storage.list_objects(pattern)
        return [f for f in all_files if f.endswith('.txt')]

    def _get_info_files(self, site: str) -> List[str]:
        """Get all info files for a site from cloud storage."""
        rx_prefix = self.config.aws_rx_prefix
        if site == 'all':
            pattern = f"{rx_prefix}"
        else:
            pattern = f"{rx_prefix}fpi"

        # Get all txt files
        all_files = self.cloud_storage.list_objects(pattern)
        txt_files = [f for f in all_files if f.endswith('.txt') and site in f]

        # Filter out cloud sensor files
        return [f for f in txt_files if not os.path.basename(f).startswith('Cloud_')]

    def _parse_info_file(self, info_file: str) -> FileMetadata:
        """Parse the info file from cloud storage to extract metadata."""
        content = self.cloud_storage.read_text_file(info_file)
        if not content:
            raise ProcessingError(f"Failed to read info file: {info_file}", {"file": info_file})
            
        lines = content.strip().split('\n')
        if len(lines) < 4:
            raise ProcessingError(f"Invalid info file format: {info_file}", {"file": info_file})
            
        filename = lines[0].rstrip().split('.tar.gz', 1)[0] + '.tar.gz'
        parts = int(lines[1])
        size = int(lines[2])
        timestamp = datetime.strptime(lines[3][:19], '%Y-%m-%d %H:%M:%S')
        
        instr_inum, site, dates = filename.split('_')
        instrument = instr_inum[0:3].lower()
        instrument_number = instr_inum[3:5]
        
        return FileMetadata(
            instrument=instrument,
            site=site,
            date=timestamp,
            instrument_number=instrument_number,
            original_filename=filename,
            expected_parts=parts,
            total_size=size,
            emails=activeinstruments()[site][instrument][instrument_number]['email']
        )

    def _acquire_lock(self, site: str):
        """Acquire processing lock for a site."""
        pid = str(os.getpid())
        pidfile = f"/tmp/Sorter_{site}.pid"
        
        if os.path.exists("/tmp/Sorterall.pid"):
            raise ProcessingError("Sorterall is running", {"site": site})
        if os.path.exists(pidfile):
            raise ProcessingError(f"Another instance is processing {site}", {"site": site})
            
        with open(pidfile, 'w') as f:
            f.write(pid)

    def _release_lock(self, site: str):
        """Release processing lock for a site."""
        pidfile = f"/tmp/Sorter_{site}.pid"
        if os.path.exists(pidfile):
            os.unlink(pidfile)

    def _handle_error(self, error: Exception, info_file: str):
        """Handle processing errors."""
        self.logger.error(f"Error processing {info_file}: {str(error)}", exc_info=True)
        
        # If this is a critical error, consider creating a GitHub issue
        if isinstance(error, ProcessingError):
            # Implement further error handling if needed
            pass

def main():
    print("Starting")


    load_dotenv()

    parser = OptionParser(usage="usage: Sorter -s SITE")
    parser.add_option("-s", "--site", dest="site", 
                     help="specific site to process (eg site,other,all)",
                     metavar="SITE", type="str", default='all')
    options, _ = parser.parse_args()
    
    print(f"Starting Sorter with site: {options.site}")
    
    try:
        logger = setup_logging(options.site)
        print(f"Logging setup complete, log directory: {'/home/jmakela/logs'}")
    except Exception as e:
        print(f"Error setting up logging: {str(e)}")
        raise
    
    config = Configuration()
    processor = DataProcessor(config)
    
    if options.site.lower() == 'all':
        sites = list(fpiinfo.get_all_sites_info().keys())
    else:
        sites = [options.site.lower()]
    
    for site in sites:
        processor.process_site(site)

if __name__ == "__main__":
    main()
