from dagster import build_asset_context

from airglow.dagster_airglow.analysis_asset import AnalysisConfig, analyze_data


def test_analysis(s3_resource):
    context = build_asset_context(
        resources={"s3": s3_resource}
    )
    config = AnalysisConfig()
    analyze_data(context, config)