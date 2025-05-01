#!/usr/bin/python
'''
Script to send files back to Airglow
  -s to add SITE

History: 27 Sep 2012 - initial script written; Daniel J. Fisher (dfisher2@illinois.edu)
         12 Feb 2014 - Updated to v3.0 - txtcheck; Daniel J. Fisher (dfisher2@illinois.edu)

'''

# Import required modules
import os
from glob import glob
import boto3
from botocore.config import Config
from botocore.exceptions import ClientError, NoCredentialsError
from dotenv import load_dotenv

def main():
    # Load environment variables from .env file
    load_dotenv()

    # Set up correct folder options
    folders = {
#        '/cygdrive/c/Sending/',
#        '/cygdrive/d/Sending/',
#        '/cygdrive/f/Sending/',
#        '/home/gps/Sending/',
#        '/home/airglow/Sending/',
#        '/home/scintmon/Sending/',
#        '/data/Sending/',
#        'C:/Sending/',
#        'D:/Sending/',
#        'F:/Sending/'
        '/home/airglow/airglow/Sending/',
        '/data/Sending/'
    }

    # destiny = 'tx@remote2.ece.illinois.edu:/rdata/airglow/rx/.'
    destiny = 'tx@remote2.ece.illinois.edu:/home/tx/rx/.'

    mfs = 70  # minimum file size to send

    send(folders, destiny, mfs)

def s3_client():
    required_vars = ['AWS_ACCESS_KEY_ID', 'AWS_SECRET_ACCESS_KEY',
                     "DEST_BUCKET", "SITE_PATH"]
    missing_vars = [var for var in required_vars if not os.environ.get(var)]

    # If any required variables are missing, log and skip upload
    if missing_vars:
        print("Skipping object store upload: Missing environment variables: ", missing_vars)
        return False

    # Check for endpoint URL (optional)
    endpoint_url = os.environ.get('AWS_S3_ENDPOINT_URL')

    # Create a boto3 client with the appropriate configuration
    try:
        s3_config = Config(
            signature_version='s3',
            s3={'addressing_style': 'path'}
        )

        s3_client_args = {
            'service_name': 's3',
            'aws_access_key_id': os.environ.get('AWS_ACCESS_KEY_ID'),
            'aws_secret_access_key': os.environ.get('AWS_SECRET_ACCESS_KEY'),
            'config': s3_config
        }

        # Add region if specified
        if os.environ.get('AWS_REGION'):
            s3_client_args['region_name'] = os.environ.get('AWS_REGION')

        # Add endpoint URL if specified
        if endpoint_url:
            s3_client_args['endpoint_url'] = endpoint_url

        return boto3.client(**s3_client_args)

    except NoCredentialsError:
        print("Credentials not available or invalid")
        return False


def upload_file_to_s3(s3, file_path: str, bucket_name: str, object_name=None):
    """
    Upload a file to an S3 bucket using environment variables for credentials

    :param file_path: File to upload
    :param bucket_name: Bucket to upload to
    :param object_name: S3 object name. If not specified, file_path's basename is used
    :return: True if file was uploaded, False otherwise
    """

    site_path = os.environ.get("SITE_PATH").rstrip("/")

    # If object_name was not specified, use file_path's basename
    if object_name is None:
        object_name = f"{site_path}/{os.path.basename(file_path)}"

    # Check if the file exists
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return False

    try:
        print(f"Starting upload of {file_path} to {bucket_name}/{object_name}")
        s3.upload_file(file_path, bucket_name, object_name)
        return True

    except ClientError as e:
        print(f"S3 client error: {e}")
        return False
    except Exception as e:
        print(f"Unexpected error during upload: {e}")
        return False


def send(folders: set[str], destiny: str, mfs: int):
    s3 = s3_client()
    bucket_name = os.environ.get("DEST_BUCKET", None)

    # Go through all possible folders for data
    files = []
    for f in folders:
        f = f+"/" if f[-1] != "/" else f
        files = files + glob(f+'*.tar.gz*')
        files = files + glob(f+'*.txt')

    # Send all files one by one; remove if sent successfully
    for f in [fx for fx in files if os.stat(fx).st_size > mfs]:
        print(f)
        if s3 and bucket_name:
            s3_result = upload_file_to_s3(s3, f, bucket_name)
        else:
            s3_result = True

#        os.system('chmod 774 '+f)
#        flag = os.system('scp ' + f + ' ' + destiny)
#        if flag == 0 and s3_result:
        if s3_result:
            os.remove(f)
            print('Completed...')

    print('All Sending is now Complete!')

if __name__ == "__main__":
    main()
