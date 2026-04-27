# download_manager.py
import os
import requests
import hashlib
import gzip
import shutil
from pathlib import Path
from typing import Optional, Dict
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
import threading
import queue
import time

class DownloadStatus(Enum):
    """下载任务状态"""
    PENDING = "pending"
    DOWNLOADING = "downloading"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"

@dataclass
class DownloadTask:
    """下载任务数据结构"""
    task_id: str
    dataset_id: str
    download_url: str
    local_path: str
    status: DownloadStatus
    progress: float  # 0-100
    file_size: Optional[int]
    downloaded_size: int
    error_message: Optional[str]
    created_at: datetime
    completed_at: Optional[datetime]

    def to_dict(self) -> Dict:
        return {
            "task_id": self.task_id,
            "dataset_id": self.dataset_id,
            "status": self.status.value,
            "progress": f"{self.progress:.1f}%",
            "file_size": f"{self.file_size / (1024**3):.2f} GB" if self.file_size else "Unknown",
            "downloaded": f"{self.downloaded_size / (1024**3):.2f} GB",
            "local_path": self.local_path if self.status == DownloadStatus.COMPLETED else None,
            "error": self.error_message,
            "created_at": self.created_at.strftime("%Y-%m-%d %H:%M:%S"),
            "completed_at": self.completed_at.strftime("%Y-%m-%d %H:%M:%S") if self.completed_at else None
        }

class DownloadManager:
    """异步下载管理器"""

    # OpenGWAS 数据下载 URL 模板
    OPENGWAS_URL_TEMPLATE = "https://gwas.mrcieu.ac.uk/files/{dataset_id}/{dataset_id}.vcf.gz"
    
    # EBI GWAS Catalog URL 模板
    EBI_URL_TEMPLATE = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{accession}"

    def __init__(self, download_dir: str = "./gwas_data"):
        self.download_dir = Path(download_dir)
        self.download_dir.mkdir(parents=True, exist_ok=True)
        
        # 任务队列和状态追踪
        self.tasks: Dict[str, DownloadTask] = {}
        self.task_queue = queue.Queue()
        self.lock = threading.Lock()
        
        # 启动后台下载线程
        self.worker_thread = threading.Thread(target=self._download_worker, daemon=True)
        self.worker_thread.start()

    def _generate_task_id(self, dataset_id: str) -> str:
        """生成唯一任务 ID"""
        timestamp = datetime.now().isoformat()
        return hashlib.md5(f"{dataset_id}_{timestamp}".encode()).hexdigest()[:12]

    def _resolve_download_url(
        self,
        dataset_id: str,
        download_source: str,
        custom_url: Optional[str] = None
    ) -> str:
        """解析下载 URL"""
        if download_source == "Custom_URL":
            if not custom_url:
                raise ValueError("Custom_URL 需要提供 custom_url 参数")
            return custom_url
        
        elif download_source == "OpenGWAS":
            return self.OPENGWAS_URL_TEMPLATE.format(dataset_id=dataset_id)
        
        elif download_source == "EBI_FTP":
            accession = dataset_id.split('-')[-1] if '-' in dataset_id else dataset_id
            return self.EBI_URL_TEMPLATE.format(accession=accession)
        
        else:
            raise ValueError(f"不支持的下载源: {download_source}")

    def submit_task(
        self,
        dataset_id: str,
        download_source: str = "OpenGWAS",
        custom_url: Optional[str] = None
    ) -> str:
        """
        提交下载任务

        Args:
            dataset_id: 数据集 ID
            download_source: 下载源
            custom_url: 自定义 URL（当 download_source 为 Custom_URL 时必填）

        Returns:
            任务 ID
        """
        try:
            # 解析下载 URL
            download_url = self._resolve_download_url(dataset_id, download_source, custom_url)
            
            # 生成任务 ID 和本地路径
            task_id = self._generate_task_id(dataset_id)
            filename = f"{dataset_id}.tsv.gz"
            local_path = str(self.download_dir / filename)

            # 创建任务对象
            task = DownloadTask(
                task_id=task_id,
                dataset_id=dataset_id,
                download_url=download_url,
                local_path=local_path,
                status=DownloadStatus.PENDING,
                progress=0.0,
                file_size=None,
                downloaded_size=0,
                error_message=None,
                created_at=datetime.now(),
                completed_at=None
            )

            # 注册任务
            with self.lock:
                self.tasks[task_id] = task
                self.task_queue.put(task_id)

            print(f"✅ 下载任务已提交: {task_id}")
            print(f"   数据集: {dataset_id}")
            print(f"   下载源: {download_source}")
            
            return task_id

        except Exception as e:
            print(f"❌ 提交下载任务失败: {str(e)}")
            raise

    def get_task_status(self, task_id: str) -> Optional[Dict]:
        """查询任务状态"""
        with self.lock:
            task = self.tasks.get(task_id)
            if not task:
                return None
            return task.to_dict()

    def _download_worker(self):
        """后台下载工作线程"""
        while True:
            try:
                task_id = self.task_queue.get(timeout=1)
                self._execute_download(task_id)
            except queue.Empty:
                continue
            except Exception as e:
                print(f"❌ 下载工作线程异常: {str(e)}")

    def _execute_download(self, task_id: str):
        """执行实际下载"""
        with self.lock:
            task = self.tasks.get(task_id)
            if not task:
                return
            task.status = DownloadStatus.DOWNLOADING

        try:
            print(f"🔄 开始下载: {task.dataset_id}")
            
            # 发起 HTTP 请求
            response = requests.get(
                task.download_url,
                stream=True,
                timeout=30,
                headers={'User-Agent': 'OpenClaw-GWAS-Agent/1.0'}
            )
            response.raise_for_status()

            # 获取文件大小
            file_size = int(response.headers.get('content-length', 0))
            with self.lock:
                task.file_size = file_size

            # 流式下载
            downloaded = 0
            chunk_size = 8192
            
            with open(task.local_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        
                        # 更新进度
                        with self.lock:
                            task.downloaded_size = downloaded
                            if file_size > 0:
                                task.progress = (downloaded / file_size) * 100

            # 如果是 gz 文件，尝试解压
            if task.local_path.endswith('.gz'):
                self._decompress_file(task)

            # 标记完成
            with self.lock:
                task.status = DownloadStatus.COMPLETED
                task.progress = 100.0
                task.completed_at = datetime.now()

            print(f"✅ 下载完成: {task.dataset_id}")
            print(f"   保存路径: {task.local_path}")

        except requests.exceptions.RequestException as e:
            with self.lock:
                task.status = DownloadStatus.FAILED
                task.error_message = f"网络错误: {str(e)}"
            print(f"❌ 下载失败: {task.dataset_id} - {str(e)}")

        except Exception as e:
            with self.lock:
                task.status = DownloadStatus.FAILED
                task.error_message = str(e)
            print(f"❌ 下载异常: {task.dataset_id} - {str(e)}")

    def _decompress_file(self, task: DownloadTask):
        """解压 gzip 文件"""
        try:
            decompressed_path = task.local_path.replace('.gz', '')
            print(f"🔄 正在解压: {task.local_path}")
            
            with gzip.open(task.local_path, 'rb') as f_in:
                with open(decompressed_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            # 删除压缩文件，更新路径
            os.remove(task.local_path)
            task.local_path = decompressed_path
            
            print(f"✅ 解压完成: {decompressed_path}")
        
        except Exception as e:
            print(f"⚠️ 解压失败（保留压缩文件）: {str(e)}")


# Function Schemas for LLM
SUBMIT_DOWNLOAD_SCHEMA = {
    "name": "submit_download_task",
    "description": "提交一个后台异步下载任务，用于下载 GWAS 汇总统计数据文件。任务提交后立即返回任务ID，下载在后台进行。",
    "parameters": {
        "type": "object",
        "properties": {
            "dataset_id": {
                "type": "string",
                "description": "需要下载的 GWAS 数据集 ID，例如 'ebi-a-GCST006867' 或 'ieu-a-2'。"
            },
            "download_source": {
                "type": "string",
                "enum": ["OpenGWAS", "EBI_FTP", "Custom_URL"],
                "description": "指定下载源。OpenGWAS 适用于 IEU 数据库，EBI_FTP 适用于 GWAS Catalog。",
                "default": "OpenGWAS"
            },
            "custom_url": {
                "type": "string",
                "description": "当 download_source 为 'Custom_URL' 时，在此填入完整的下载链接。"
            }
        },
        "required": ["dataset_id", "download_source"]
    }
}

CHECK_STATUS_SCHEMA = {
    "name": "check_download_status",
    "description": "查询后台下载任务的当前状态，包括进度、文件大小、是否完成等信息。",
    "parameters": {
        "type": "object",
        "properties": {
            "task_id": {
                "type": "string",
                "description": "由 submit_download_task 返回的任务 ID。"
            }
        },
        "required": ["task_id"]
    }
}


# 测试代码
if __name__ == "__main__":
    manager = DownloadManager(download_dir="./test_downloads")

    # 测试提交任务
    print("=== 测试提交下载任务 ===")
    task_id = manager.submit_task(
        dataset_id="ieu-a-2",
        download_source="OpenGWAS"
    )

    # 轮询任务状态
    print("\n=== 监控下载进度 ===")
    while True:
        status = manager.get_task_status(task_id)
        if status:
            print(f"状态: {status['status']} | 进度: {status['progress']}")
            
            if status['status'] in ['completed', 'failed']:
                print("\n最终状态:")
                for key, value in status.items():
                    print(f"  {key}: {value}")
                break
        
        time.sleep(2)
