# Airglow Dagster Pipeline

This Dagster project manages the processing pipeline for Airglow instrument data. It automatically processes new observations uploaded to an S3 bucket, performs analysis, and cleans up raw data files.

## Components

### Sensor
- `instrument_upload_sensor`: Monitors an S3 bucket for new instrument data uploads
  - Runs every hour
  - Detects new chunked archive files and associated cloud cover data
  - Triggers the analysis job when complete data sets are available

### Assets
1. `unzip_chunked_archive`
   - Processes chunked tar.gz archives uploaded to S3
   - Combines chunks into a single archive
   - Extracts files to the appropriate S3 location
   - Copies cloud cover data to the archive directory
   - Returns configuration for subsequent analysis steps

2. `analyze_data_pipeline`
   - Downloads FPI (Fabry-Perot Interferometer) data and cloud cover data
   - Performs analysis on the instrument data
   - Uploads analysis results back to S3
   - Stores results in MySQL database

3. `delete_raw`
   - Cleans up raw data files after successful processing
   - Removes chunked archives, cloud cover files, and instrument log files
   - Only executes after successful analysis

### Job
- `analysis_job`: Orchestrates the complete data processing workflow
  - Executes assets in sequence: unzip → analyze → delete
  - Triggered by the instrument upload sensor
  - Processes data for specific sites and observation dates

## Data Flow
1. New instrument data is uploaded to S3 as chunked archives
2. The sensor detects complete data sets (archive chunks + cloud cover data)
3. The analysis job is triggered with appropriate configuration
4. Raw data is processed and analyzed
5. Results are stored in S3 and MySQL
6. Raw data files are cleaned up

## Configuration
The pipeline requires the following environment variables:
- `DEST_BUCKET`: S3 bucket name for data storage
- MySQL database credentials (configured through Dagster resources)

## Dependencies
- dagster
- dagster-aws
- dagster-mysql
- dagster-ncsa
- airglow (internal package)

## Support

### Reanalyzing Data
The pipeline includes a `reanalyze_data` asset that allows you to re-run the analysis for a specific instrument and date without needing to re-extract the data. This is useful for:
- Fixing analysis issues
- Updating results with new analysis parameters
- Recovering from failed analysis runs

To reanalyze data:

1. Navigate to the Dagster UI at https://airglow.software.ncsa.illinois.edu/assets/reanalyze_data
2. Click the "Materialize" dropdown button
3. Select "Open Launchpad"
4. In the Launchpad, provide the following configuration:
   ```yaml
   config:
     site: "uao"  # Site code (e.g., uao)
     year: "2025"  # Year of observation
     observation_date: "20250429"  # Date in YYYYMMDD format
     fpi_data_path: "fpi/minime05/uao/2025/20250429"  # Path to FPI data in S3
     cloud_cover_path: "cloudsensor/uao/2025"  # Path to cloud cover data in S3
   ```
5. Click "Launch Run" to start the reanalysis

The reanalysis will use the existing data in S3 and MySQL, perform the analysis again, and update the results. 