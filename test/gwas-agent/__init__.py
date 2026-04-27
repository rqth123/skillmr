# gwas_agent/__init__.py
"""
GWAS Analysis Agent for OpenClaw
专注于 GWAS 数据分析的 AI Agent
"""

__version__ = "1.0.0"
__author__ = "LYJJ"

from .agent import GWASAgent
from .skills import (
    search_gwas_catalog,
    submit_download_task,
    check_download_status,
    peek_data_headers,
    validate_gwas_format,
    execute_analysis_code
)

__all__ = [
    'GWASAgent',
    'search_gwas_catalog',
    'submit_download_task',
    'check_download_status',
    'peek_data_headers',
    'validate_gwas_format',
    'execute_analysis_code'
]
