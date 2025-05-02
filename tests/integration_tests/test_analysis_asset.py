import logging
import sys
from dagster import build_asset_context

from airglow.dagster_airglow.assets.analysis_asset import AnalysisConfig, reanalyze_data

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    stream=sys.stdout,  # Explicitly set to stdout
    force=True,  # Force configuration even if already configured
)
logger = logging.getLogger(__name__)


def test_analysis(s3_resource):
    context = build_asset_context(resources={"s3": s3_resource, "mysql": None})
    config = AnalysisConfig(
        observation_date="20250425",
        year="2025",
        site="uao",
        fpi_data_path="fpi/minime05/uao/2025/20250425",
        cloud_cover_path="cloudsensor/uao",
    )
    reanalyze_data(
        context,
        config,
    )
