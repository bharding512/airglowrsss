#!/usr/bin/env python3
"""
MANGO Data Pipeline Sorter

This module handles the sorting and processing of instrument data files,
primarily focused on FPI and Cloud sensor data with extensibility for other instruments.
"""

import os
import sys
import shutil
import tarfile
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Union
from dataclasses import dataclass
import logging
from abc import ABC, abstractmethod
import numpy as np
import pytz
from glob import glob

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

def setup_logging(site: str):
    """Configure logging for the Sorter."""
    # Create logs directory if it doesn't exist
    log_dir = '/home/airglow/logs'
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

class Config:
    """Configuration management."""
    def __init__(self):
        self.rx_dir = '/home/tx/rx/'
        self.data_dir = '/rdata/airglow/'
        self.share_dir = '/rdata/airglow/share/'
        
        # Processing configurations
        self.fpi_process_kwargs = {
            'reference': 'laser',
            'send_to_website': True,
            'enable_share': False,
            'send_to_madrigal': True,
            'sky_line_tag': 'X',
            'fpi_dir': '/rdata/airglow/fpi/',
            'bw_dir': '/rdata/airglow/templogs/cloudsensor/',
            'x300_dir': '/rdata/airglow/templogs/x300/',
            'results_stub': '/rdata/airglow/fpi/results/'
        }

class InstrumentHandler(ABC):
    """Base class for instrument-specific handlers."""
    def __init__(self, config: Config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
    @abstractmethod
    def sort_files(self, metadata: FileMetadata) -> List[str]:
        """Sort files for this instrument type."""
        pass
        
    @abstractmethod
    def process_data(self, metadata: FileMetadata, sorted_files: List[str]):
        """Process the sorted data."""
        pass

    def get_target_directory(self, metadata: FileMetadata) -> Path:
        """Get the target directory for the instrument data."""
        return Path(self.config.data_dir)

class FPIHandler(InstrumentHandler):
    """Handler for FPI instruments."""
    def __init__(self, config: Config):
        super().__init__(config)
        self.issue_manager = SiteIssueManager(
            token=github_config.github_token,
            repo_name=github_config.github_repo
        )

    def get_target_directory(self, metadata: FileMetadata) -> Path:
        """Get the target directory for FPI data."""
        return Path(self.config.data_dir) / 'fpi' / f"minime{metadata.instrument_number}" / metadata.site / str(metadata.date.year)

    def sort_files(self, metadata: FileMetadata) -> List[str]:
        target_dir = self.get_target_directory(metadata)
        self.logger.info(f"Sorting files to: {target_dir}")
        
        target_dir.mkdir(parents=True, exist_ok=True)
        result = self._concat_and_extract(metadata, target_dir)
        
        if result:
            self.logger.info(f"Successfully sorted {len(result)} files")
            self._set_permissions(target_dir, result)
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
        finally:
            self.current_metadata = None

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
        """Concatenate split files and extract contents."""
        file_pattern = f"{metadata.original_filename}*"
        parts = glob(str(Path(self.config.rx_dir) / file_pattern))
        
        if len(parts) != metadata.expected_parts:
            self.logger.warning(f"Expected {metadata.expected_parts} parts, found {len(parts)}")
            return []

        temp_file = target_dir / f"temp_{datetime.now().strftime('%Y%m%d%H%M%S%f')}.tar.gz"
        try:
            # Concatenate parts
            with open(temp_file, 'wb') as outfile:
                for part in sorted(parts):
                    with open(part, 'rb') as infile:
                        outfile.write(infile.read())

            # Verify size
            if temp_file.stat().st_size != metadata.total_size:
                self.logger.error("Size mismatch after concatenation")
                return []

            # Extract
            with tarfile.open(temp_file, 'r:gz') as tar:
                tar.extractall(path=target_dir)
                extracted_files = [member.name for member in tar.getmembers()]

            # Clean up parts after successful extraction
            for part in parts:
                os.remove(part)
            self.logger.info(f"Cleaned up {len(parts)} partial files")

            return extracted_files

        finally:
            if temp_file.exists():
                temp_file.unlink()

    def _get_local_datetime(self, first_file: str) -> datetime:
        """Get local datetime from first file."""
        full_path = self.get_target_directory(self.current_metadata) / first_file
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
            shutil.chown(file_path, 'airglow', 'fpi')

    def _handle_processing_error(self, error: Exception, metadata: FileMetadata, 
                               year: int, doy: int, tag: str):
        """Handle processing errors."""
        self.issue_manager.handle_processing_issue(
            site_id=metadata.site,
            message=str(error),
            category=IssueType.ERROR,
            error_type=type(error).__name__,
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
    def sort_files(self, metadata: FileMetadata) -> List[str]:
        target_dir = self.get_target_directory(metadata) / 'templogs/cloudsensor' / metadata.site
        target_dir.mkdir(parents=True, exist_ok=True)
        
        source_file = Path(self.config.rx_dir) / metadata.original_filename
        target_file = target_dir / metadata.original_filename
        shutil.copy2(source_file, target_file)
        
        return [str(target_file)]

    def process_data(self, metadata: FileMetadata, sorted_files: List[str]):
        """Process cloud sensor data - implement specific processing if needed."""
        pass

    """Handler for cloud sensor data."""
    def get_target_directory(self, site: str) -> Path:
        """Get the target directory for cloud sensor data."""
        return Path(self.config.data_dir) / 'templogs/cloudsensor' / site

    def process_cloud_file(self, cloud_file: Path):
        """Process a single cloud sensor file."""
        # Extract site from filename (Cloud_site_date.txt)
        _, site, _ = cloud_file.stem.split('_')
        
        # Create target directory
        target_dir = self.get_target_directory(site)
        target_dir.mkdir(parents=True, exist_ok=True)
        
        # Move file to target directory, overwriting if exists
        target_file = target_dir / cloud_file.name
        shutil.move(str(cloud_file), str(target_file))
        os.chmod(target_file, 0o644)
        self.logger.info(f"Moved cloud sensor file to {target_file}")

class DataProcessor:
    """Main data processing orchestrator."""
    def __init__(self, config: Config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.handlers = {
            'fpi': FPIHandler(config),
            'cloud': CloudSensorHandler(config),
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
                    handler = self.handlers.get(metadata.instrument)
                    if handler:
                        self.logger.info(f"Using handler for instrument: {metadata.instrument}")
                        sorted_files = handler.sort_files(metadata)
                        if sorted_files:
                            self.logger.info(f"Sorted {len(sorted_files)} files")
                            os.remove(info_file)
                            handler.process_data(metadata, sorted_files)
                            self.logger.info("Processing complete")
                    else:
                        self.logger.warning(f"No handler for instrument: {metadata.instrument}")
                except Exception as e:
                    self.logger.error(f"Error processing {info_file}: {str(e)}", exc_info=True)
                    self._handle_error(e, info_file)
        finally:
            self._release_lock(site)
            self.logger.info(f"Completed processing for site: {site}")

    def get_cloud_files(self, site: str) -> List[Path]:
        """Get cloud sensor files for a site."""
        rx_path = Path(self.config.rx_dir)
        if site == 'all':
            pattern = 'Cloud_*.txt'
        else:
            pattern = f'Cloud_{site}_*.txt'
        return list(rx_path.glob(pattern))

    def _get_info_files(self, site: str) -> List[Path]:
        """Get all info files for a site."""
        rx_path = Path(self.config.rx_dir)
        if site == 'all':
            pattern = '*.txt'
        else:
            pattern = f'*_{site}_*.txt'

        # Get all txt files
        all_files = list(rx_path.glob(pattern))

        # Filter out cloud sensor files
        return [f for f in all_files if not f.name.startswith('Cloud_')]

    def _parse_info_file(self, info_file: Path) -> FileMetadata:
        """Parse the info file to extract metadata."""
        with open(info_file, 'r') as f:
            lines = f.readlines()
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

    def _handle_error(self, error: Exception, info_file: Path):
        """Handle processing errors."""
        self.logger.error(f"Error processing {info_file}: {str(error)}", exc_info=True)

def main():
    print("Starting")
    parser = OptionParser(usage="usage: Sorter -s SITE")
    parser.add_option("-s", "--site", dest="site", 
                     help="specific site to process (eg site,other,all)",
                     metavar="SITE", type="str", default='all')
    options, _ = parser.parse_args()
    
    print(f"Starting Sorter with site: {options.site}")
    
    try:
        logger = setup_logging(options.site)
        print(f"Logging setup complete, log directory: {'/home/airglow/logs'}")
    except Exception as e:
        print(f"Error setting up logging: {str(e)}")
        raise
    
    config = Config()
    processor = DataProcessor(config)
    
    if options.site.lower() == 'all':
        sites = list(fpiinfo.get_all_sites_info().keys())
    else:
        sites = [options.site.lower()]
    
    for site in sites:
        processor.process_site(site)

if __name__ == "__main__":
    main()
