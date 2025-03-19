# AirGlow-RSSS
Tools to support the understanding the dynamics of the ionosphere using a
combination of optical and radio techniques. This is accomplished through 
instrument development, experimental campaigns, and data analysis. 
[Our group](https://airglow.ece.illinois.edu/) is 
part of the Remote Sensing & Space Science group in the Department of Electrical 
and Computer Engineering. We are always looking for students interested in becoming 
involved in our research, especially those with a background in signal/image 
processing, plasma physics, or optical/radio instrument design.

## Sender
This script is used to send data files to the servers at Illinois. As part of a migration
away from custom servers, the script will optionally upload the files to the projects Open 
Storage Network (OSN) bucket.

### Installation
The Sender script is a Python script that requires Python 3.10 or higher. It can 
be installed using pip:

```bash
 pip install git+https://github.com/AirglowRSSS/airglowrsss.git#subdirectory=DataManagement
```

### Configuration
In order to upload files to S3, the script requires AWS credentials. You can set these
as environment variables:

| Environment Variable  | Description                                     | Example Value                              |
|-----------------------|-------------------------------------------------|--------------------------------------------|
| AWS_ACCESS_KEY_ID     | Access key for AWS authentication               | `AKIAIOSFODNN7EXAMPLE`                     |
| AWS_SECRET_ACCESS_KEY | Secret key for AWS authentication               | `wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY` |
| AWS_S3_ENDPOINT_URL   | Custom endpoint URL for S3 API                  | `https://ncsa.osn.xsede.org`               |
| AWS_S3_REGION         | AWS region for S3 bucket  (only needed for AWS) | `us-east-1`                                |
| DEST_BUCKET           | Destination S3 bucket name                      | `my-backup-bucket`                         |
| SITE_PATH             | Prfix to be applied to uploaded object name     | `airglow/rar/oau`                          |

These can be set in the shell, or can be put in a `.env` file in the same directory as the script.

### Usage
To use the script, you can run it from the command line with this simple command:

```bash
% Sender
```

