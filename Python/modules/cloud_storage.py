#!/usr/bin/env python3
# cloud_storage.py
# Classes for AWS cloud storage operations and configuration management

import boto3
from botocore.exceptions import ClientError
from botocore.config import Config
from dotenv import load_dotenv
from typing import Dict, List, Optional, Union, Tuple
import io
import os
import logging


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
        """List objects in the bucket with the given prefix, excluding the directory itself."""
        try:
            response = self.client_s3v4.list_objects_v2(
                Bucket=self.bucket,
                Prefix=prefix
            )
            
            if 'Contents' in response:
                # Filter out objects that match the prefix exactly (which is the directory itself)
                return [obj['Key'] for obj in response['Contents'] if obj['Key'] != prefix]
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
    """Configuration management for AWS cloud storage and FPI processing."""
    def __init__(self, env_file='.env'):
        # Load environment variables from .env file
        load_dotenv(env_file)
        
        # AWS Cloud Configuration
        self.aws_access_key = os.getenv('AWS_ACCESS_KEY')
        self.aws_secret_key = os.getenv('AWS_SECRET_KEY')
        self.aws_endpoint_url = os.getenv('AWS_ENDPOINT_URL', 'https://ncsa.osn.xsede.org')
        self.dest_bucket = os.getenv('DEST_BUCKET', 'airglow')
        
        # Cloud directory structure
        self.aws_rx_prefix = os.getenv('AWS_RX_PREFIX', 'raw/')
        self.aws_data_prefix = os.getenv('AWS_DATA_PREFIX', '')
        self.aws_results_prefix = os.getenv('AWS_RESULTS_PREFIX', 'results/')
        self.aws_summaryimages_prefix = os.getenv('AWS_SUMMARYIMAGES_PREFIX', 'SummaryImages/')
        self.aws_madrigal_prefix = os.getenv('AWS_MADRIGAL_PREFIX', 'madrigal/')
        self.aws_fpi_prefix = os.getenv('AWS_FPI_PREFIX', 'fpi/')
        self.aws_cloudsensor_prefix = os.getenv('AWS_CLOUDSENSOR_PREFIX', 'cloudsensor/')
        self.aws_mango_movies_prefix = os.getenv('AWS_MANGO_MOVIES_PREFIX', 'mango_movies/')
        
        # Local temporary directory
        self.temp_dir = os.getenv('TEMP_DIR', '/home/jmakela/tmp/mango')
        os.makedirs(self.temp_dir, exist_ok=True)
        
        # Processing configurations
        self.fpi_process_kwargs = {
            'reference': os.getenv('FPI_REFERENCE', 'laser'),
            'send_to_website': os.getenv('FPI_SEND_TO_WEBSITE', 'True').lower() == 'true',
            'enable_share': os.getenv('FPI_ENABLE_SHARE', 'False').lower() == 'true',
            'send_to_madrigal': os.getenv('FPI_SEND_TO_MADRIGAL', 'True').lower() == 'true',
            'sky_line_tag': os.getenv('FPI_SKY_LINE_TAG', 'XR'),
            'fpi_dir': os.getenv('FPI_DIR', f"{self.temp_dir}/fpi/"),
            'bw_dir': os.getenv('BW_DIR', f"{self.temp_dir}/cloudsensor/"),
            'x300_dir': os.getenv('X300_DIR', f"{self.temp_dir}/templogs/x300/"),
            'results_stub': os.getenv('RESULTS_STUB', f"{self.temp_dir}/results/"),
            'temp_plots_stub': os.getenv('TEMP_PLOTS_STUB', f"{self.temp_dir}/temporary_plots/"),
            'madrigal_stub': os.getenv('MADRIGAL_STUB', f"{self.temp_dir}/madrigal/"),
        }

        # Create necessary directories
        for directory in [
            self.fpi_process_kwargs['fpi_dir'],
            self.fpi_process_kwargs['bw_dir'],
            self.fpi_process_kwargs['x300_dir'],
            self.fpi_process_kwargs['results_stub'],
            self.fpi_process_kwargs['temp_plots_stub'],
            self.fpi_process_kwargs['madrigal_stub']
        ]:
            os.makedirs(directory, exist_ok=True)
