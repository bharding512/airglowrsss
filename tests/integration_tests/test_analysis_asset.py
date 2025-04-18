from dagster import build_asset_context

from airglow.dagster_airglow.analysis_asset import analyze_data


def test_analysis(s3_resource):
    context = build_asset_context(
        resources={"s3": s3_resource}
    )
    analyze_data(context,
       {
            "observation_date": "20250416",
            "year": "2025",
            "site": "blo",
            "fpi_data_path": "fpi/minime12/blo/2025/20250416",
            "cloud_cover_path": "cloudsensor/blo/2025"},
                 sky_line_tag="XR", s3=s3_resource)