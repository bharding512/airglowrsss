from unittest.mock import MagicMock

from dagster_airglow.sensors import (cloud_cover_files_for_site,
                                     group_files_by_date,
                                     instrument_upload_sensor)
import dagster as dg


def test_group_files_by_date():
    # Test data
    files = [
        "raw/Cloud_ABC_20230101.txt",
        "raw/Cloud_DEF_20230101.txt",
        "raw/Cloud_XYZ_20230101.txt",
        "raw/Cloud_ABC_20230201.txt",
        "raw/Cloud_DEF_20230201.txt",
        "raw/Cloud_XYZ_20230201.txt",
        "raw/fpi05_ABC_20250409.tar.gz000006"
    ]

    # Expected output
    expected_output = {
        20230101: [
            "raw/Cloud_ABC_20230101.txt",
            "raw/Cloud_DEF_20230101.txt",
            "raw/Cloud_XYZ_20230101.txt"
        ],
        20230201: [
            "raw/Cloud_ABC_20230201.txt",
            "raw/Cloud_DEF_20230201.txt",
            "raw/Cloud_XYZ_20230201.txt"
        ],
        20250409: [
            "raw/fpi05_ABC_20250409.tar.gz000006"
        ]
    }

    # Call the function
    result = group_files_by_date(files)

    # Assert the result
    assert result == expected_output


def test_cloud_cover_files_for_site():
    # Test data
    files = [
        "raw/Cloud_ABC_20230101.txt",
        "raw/Cloud_DEF_20230101.txt",
        "raw/Cloud_XYZ_20230101.txt",
        "raw/Cloud_ABC_20230201.txt",
        "raw/Cloud_DEF_20230201.txt",
        "raw/Cloud_XYZ_20230201.txt",
        "raw/fpi05_ABC_20250409.tar.gz000006"
    ]
    site = "ABC"

    # Expected output
    expected_output = [
        "raw/Cloud_ABC_20230101.txt",
        "raw/Cloud_ABC_20230201.txt"
    ]

    # Call the function
    result = cloud_cover_files_for_site(site, files)

    # Assert the result
    assert result == expected_output


def test_instrument_upload_sensor():
    mock_s3 = MagicMock()
    paginator = MagicMock()
    files = ["raw/fpi05_ABC_20250409.tar.gz000000",
             "raw/fpi05_ABC_20250409.tar.gz000001",
             "raw/fpi05_ABC_20250409.tar.gz000002",
             "raw/fpi05_ABC_20250409.txt",
             "raw/Cloud_ABC_20250409.txt",
             "raw/Cloud_ABC_20250410.txt",
             # Incomplete upload
             "raw/fpi05_ABC_20250410.tar.gz000000",
             "raw/fpi05_ABC_20250410.tar.gz000002",
             # Just the log file
             "raw/fpi05_ABC_20250409.txt"
             ]
    paginator.paginate.return_value = [
        {"Contents": [{"Key": file} for file in files]}
    ]

    mock_s3.get_client.return_value.get_paginator.return_value = paginator

    context = dg.build_sensor_context(
        resources={"s3": mock_s3}
    )
    runs = list(instrument_upload_sensor(context))

    assert len(runs) == 1
    run09 = runs[0]
    assert run09.run_key == 'sort-20250409-ABC'
    assert run09.run_config["ops"]['unzip_chunked_archive']["config"] == {
        'observation_date': '20250409', 'site': 'ABC',
        'file_chunks': ['raw/fpi05_ABC_20250409.tar.gz000000',
                        'raw/fpi05_ABC_20250409.tar.gz000001',
                        'raw/fpi05_ABC_20250409.tar.gz000002'],
        'cloud_files': ['raw/Cloud_ABC_20250409.txt', 'raw/Cloud_ABC_20250410.txt']
    }
