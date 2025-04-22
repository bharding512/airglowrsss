from dagster import build_asset_context

from airglow.dagster_airglow.analysis_asset import analyze_data


def test_analysis(s3_resource):
    context = build_asset_context(
        resources={"s3": s3_resource}
    )
    analyze_data(context,
       {
            "observation_date": "20250420",
            "year": "2025",
            "site": "mor",
            "fpi_data_path": "fpi/minime03/mor/2025/20250420",
            "cloud_cover_path": "cloudsensor/mor"},
                 sky_line_tag="X", s3=s3_resource, mysql=None)