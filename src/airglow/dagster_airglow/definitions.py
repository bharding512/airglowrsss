
from dagster import Definitions, EnvVar, define_asset_job
from dagster_mysql import MySQLResource

from dagster_ncsa import S3ResourceNCSA

from airglow.dagster_airglow.assets.analysis_asset import analyze_data_xg, \
    analyze_data_xr, analyze_data_x
from airglow.dagster_airglow.delete_raw import delete_raw
from airglow.dagster_airglow.assets.upload_chunked_archive import unzip_chunked_archive
from airglow.dagster_airglow.sensors import instrument_upload_sensor

all_assets = [
    unzip_chunked_archive,
    analyze_data_xg,
    analyze_data_xr,
    analyze_data_x,
    delete_raw
]

analysis_job = define_asset_job(
    name="analysis_job",
    selection=[
        "unzip_chunked_archive",
        "analyze_data_xg",
        "analyze_data_xr",
        "analyze_data_x",
        "delete_raw",
    ],
)
defs = Definitions(
    assets=all_assets,
    jobs=[analysis_job],
    sensors=[instrument_upload_sensor],
    resources={
        "s3": S3ResourceNCSA(
            endpoint_url=EnvVar("AWS_S3_ENDPOINT_URL"),
            aws_access_key_id=EnvVar("AWS_ACCESS_KEY_ID"),
            aws_secret_access_key=EnvVar("AWS_SECRET_ACCESS_KEY"),
        ),
        "mysql": MySQLResource(
            host=EnvVar("MYSQL_HOST"),
            user=EnvVar("MYSQL_USER"),
            password=EnvVar("MYSQL_PASSWORD"),
            database=EnvVar("MYSQL_DATABASE"),
        )
    },
)
