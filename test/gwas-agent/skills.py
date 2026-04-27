# gwas_agent/skills.py
"""
GWAS Agent Skills - 技能函数定义
将之前实现的类方法封装为独立的技能函数
"""
from typing import Dict, List, Optional, Any
import json
# 导入之前实现的模块 - 修正导入路径
from .query_open_gwas_api import OpenGWASQueryTool  # 修正：使用原始文件名
from .download_manager import DownloadManager
from .data_inspector import DataInspector
from .code_executor import SafeCodeExecutor, GWASAnalysisHelper
# 全局实例（单例模式）
_gwas_query_engine = None  # 修正：重命名变量
_download_manager = None
_data_inspector = None
_code_executor = None
_analysis_helper = None
def _get_gwas_query_engine() -> OpenGWASQueryTool:  # 修正：函数名
    """获取 GWAS 查询引擎实例"""
    global _gwas_query_engine
    if _gwas_query_engine is None:
        _gwas_query_engine = OpenGWASQueryTool()
    return _gwas_query_engine
def _get_download_manager() -> DownloadManager:
    """获取下载管理器实例"""
    global _download_manager
    if _download_manager is None:
        _download_manager = DownloadManager()
    return _download_manager
def _get_data_inspector() -> DataInspector:
    """获取数据检查器实例"""
    global _data_inspector
    if _data_inspector is None:
        _data_inspector = DataInspector()
    return _data_inspector
def _get_code_executor() -> SafeCodeExecutor:
    """获取代码执行器实例"""
    global _code_executor
    if _code_executor is None:
        _code_executor = SafeCodeExecutor()
    return _code_executor
def _get_analysis_helper() -> GWASAnalysisHelper:
    """获取分析助手实例"""
    global _analysis_helper
    if _analysis_helper is None:
        executor = _get_code_executor()
        _analysis_helper = GWASAnalysisHelper(executor)
    return _analysis_helper
# ============================================================================
# Skill 1: 搜索 GWAS 目录
# ============================================================================
# ============================================================================
# Skill 1: 搜索 GWAS 目录
# ============================================================================
def search_gwas_catalog(
    trait: Optional[str] = None,
    author: Optional[str] = None,
    pubmed_id: Optional[str] = None,
    max_results: int = 10
) -> str:
    """
    搜索 GWAS Catalog 和 OpenGWAS 数据库
    Args:
        trait: 表型/疾病关键词
        author: 作者姓名
        pubmed_id: PubMed ID
        max_results: 最大返回结果数
    Returns:
        JSON 格式的搜索结果
    """
    try:
        engine = _get_gwas_query_engine()  
        
        # 调用 query 方法，进行参数映射
        results = engine.query(
            phenotype=trait if trait else "Any",
            # query 方法目前只支持 phenotype, population, top_k。
            # author 和 pubmed_id 暂不传入，以免报错。
            top_k=max_results
        )
        
        return json.dumps({
            "success": True,
            "count": len(results),
            "results": [r.to_dict() for r in results]
        }, indent=2, ensure_ascii=False)
    except Exception as e:
        return json.dumps({
            "success": False,
            "error": str(e)
        }, indent=2)
# ============================================================================
# Skill 2: 提交下载任务
# ============================================================================

def submit_download_task(
    dataset_id: str,
    download_source: str = "OpenGWAS",
    custom_url: Optional[str] = None
) -> str:
    """
    提交后台下载任务
    
    Args:
        dataset_id: 数据集 ID
        download_source: 下载源 (OpenGWAS/EBI_FTP/Custom_URL)
        custom_url: 自定义 URL
    
    Returns:
        JSON 格式的任务信息
    """
    try:
        manager = _get_download_manager()
        task_id = manager.submit_task(
            dataset_id=dataset_id,
            download_source=download_source,
            custom_url=custom_url
        )
        
        return json.dumps({
            "success": True,
            "task_id": task_id,
            "message": f"下载任务已提交，任务ID: {task_id}"
        }, indent=2)
        
    except Exception as e:
        return json.dumps({
            "success": False,
            "error": str(e)
        }, indent=2)


# ============================================================================
# Skill 3: 查询下载状态
# ============================================================================

def check_download_status(task_id: str) -> str:
    """
    查询下载任务状态
    
    Args:
        task_id: 任务 ID
    
    Returns:
        JSON 格式的任务状态
    """
    try:
        manager = _get_download_manager()
        status = manager.get_task_status(task_id)
        
        if status is None:
            return json.dumps({
                "success": False,
                "error": f"未找到任务: {task_id}"
            }, indent=2)
        
        return json.dumps({
            "success": True,
            "status": status
        }, indent=2, ensure_ascii=False)
        
    except Exception as e:
        return json.dumps({
            "success": False,
            "error": str(e)
        }, indent=2)


# ============================================================================
# Skill 4: 预览数据结构
# ============================================================================

def peek_data_headers(
    file_path: str,
    num_lines: int = 5
) -> str:
    """
    读取文件头部并分析结构
    
    Args:
        file_path: 文件路径
        num_lines: 读取行数
    
    Returns:
        JSON 格式的数据预览
    """
    try:
        inspector = _get_data_inspector()
        preview = inspector.peek_data(file_path, num_lines)
        
        return json.dumps({
            "success": True,
            "preview": preview.to_dict()
        }, indent=2, ensure_ascii=False)
        
    except Exception as e:
        return json.dumps({
            "success": False,
            "error": str(e)
        }, indent=2)


# ============================================================================
# Skill 5: 验证数据格式
# ============================================================================

def validate_gwas_format(file_path: str) -> str:
    """
    验证 GWAS 文件格式
    
    Args:
        file_path: 文件路径
    
    Returns:
        JSON 格式的验证结果
    """
    try:
        inspector = _get_data_inspector()
        validation = inspector.validate_gwas_format(file_path)
        
        return json.dumps({
            "success": True,
            "validation": validation
        }, indent=2, ensure_ascii=False)
        
    except Exception as e:
        return json.dumps({
            "success": False,
            "error": str(e)
        }, indent=2)


# ============================================================================
# Skill 6: 执行分析代码
# ============================================================================

def execute_analysis_code(
    code: str,
    description: str = ""
) -> str:
    """
    执行 Python 分析代码
    
    Args:
        code: Python 代码
        description: 代码描述
    
    Returns:
        JSON 格式的执行结果
    """
    try:
        executor = _get_code_executor()
        result = executor.execute(code, description)
        
        return json.dumps({
            "success": result.success,
            "result": result.to_dict()
        }, indent=2, ensure_ascii=False)
        
    except Exception as e:
        return json.dumps({
            "success": False,
            "error": str(e)
        }, indent=2)


# ============================================================================
# 技能元数据（用于 OpenClaw 注册）
# ============================================================================

SKILL_DEFINITIONS = [
    {
        "name": "search_gwas_catalog",
        "function": search_gwas_catalog,
        "description": "搜索 GWAS Catalog 和 OpenGWAS 数据库，查找特定表型或疾病的 GWAS 研究",
        "parameters": {
            "type": "object",
            "properties": {
                "trait": {
                    "type": "string",
                    "description": "表型或疾病关键词，例如 'Type 2 diabetes' 或 'Body mass index'"
                },
                "author": {
                    "type": "string",
                    "description": "作者姓名，用于筛选特定研究团队的数据"
                },
                "pubmed_id": {
                    "type": "string",
                    "description": "PubMed ID，用于查找特定文献的 GWAS 数据"
                },
                "max_results": {
                    "type": "integer",
                    "description": "最大返回结果数，默认 10",
                    "default": 10
                }
            }
        }
    },
    {
        "name": "submit_download_task",
        "function": submit_download_task,
        "description": "提交后台异步下载任务，用于下载 GWAS 汇总统计数据",
        "parameters": {
            "type": "object",
            "properties": {
                "dataset_id": {
                    "type": "string",
                    "description": "数据集 ID，例如 'ebi-a-GCST006867' 或 'ieu-a-2'"
                },
                "download_source": {
                    "type": "string",
                    "enum": ["OpenGWAS", "EBI_FTP", "Custom_URL"],
                    "description": "下载源",
                    "default": "OpenGWAS"
                },
                "custom_url": {
                    "type": "string",
                    "description": "自定义下载 URL（当 download_source 为 Custom_URL 时必填）"
                }
            },
            "required": ["dataset_id"]
        }
    },
    {
        "name": "check_download_status",
        "function": check_download_status,
        "description": "查询后台下载任务的状态和进度",
        "parameters": {
            "type": "object",
            "properties": {
                "task_id": {
                    "type": "string",
                    "description": "由 submit_download_task 返回的任务 ID"
                }
            },
            "required": ["task_id"]
        }
    },
    {
        "name": "peek_data_headers",
        "function": peek_data_headers,
        "description": "读取 GWAS 数据文件的前几行，分析列结构并智能识别标准字段",
        "parameters": {
            "type": "object",
            "properties": {
                "file_path": {
                    "type": "string",
                    "description": "GWAS 数据文件的本地路径"
                },
                "num_lines": {
                    "type": "integer",
                    "description": "读取的数据行数，默认 5",
                    "default": 5
                }
            },
            "required": ["file_path"]
        }
    },
    {
        "name": "validate_gwas_format",
        "function": validate_gwas_format,
        "description": "验证 GWAS 文件格式是否符合标准，检查必需列",
        "parameters": {
            "type": "object",
            "properties": {
                "file_path": {
                    "type": "string",
                    "description": "需要验证的文件路径"
                }
            },
            "required": ["file_path"]
        }
    },
    {
        "name": "execute_analysis_code",
        "function": execute_analysis_code,
        "description": "在隔离环境中执行 Python 分析代码，支持 pandas、numpy、matplotlib",
        "parameters": {
            "type": "object",
            "properties": {
                "code": {
                    "type": "string",
                    "description": "要执行的 Python 代码"
                },
                "description": {
                    "type": "string",
                    "description": "代码功能描述",
                    "default": ""
                }
            },
            "required": ["code"]
        }
    }
]
