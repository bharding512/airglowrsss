from dagster_airglow.sensors import cloud_cover_files_for_site


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
