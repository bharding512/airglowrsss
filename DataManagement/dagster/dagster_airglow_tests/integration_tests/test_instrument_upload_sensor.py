from dagster_airglow.sensors import instrument_upload_sensor
import dagster as dg


def test_read_raw_dataset(s3_resource):
    context = dg.build_sensor_context(
        resources={"s3": s3_resource},
        cursor="20250325"
    )
    for run in instrument_upload_sensor(context):
        print(run)
