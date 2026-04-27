# data_inspector.py
import gzip
import csv
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
from enum import Enum
import re

class FileFormat(Enum):
    """文件格式类型"""
    TSV = "tsv"
    CSV = "csv"
    VCF = "vcf"
    UNKNOWN = "unknown"

@dataclass
class ColumnMapping:
    """列名映射结果"""
    snp_col: Optional[str] = None
    chr_col: Optional[str] = None
    bp_col: Optional[str] = None
    p_col: Optional[str] = None
    beta_col: Optional[str] = None
    se_col: Optional[str] = None
    ea_col: Optional[str] = None  # Effect Allele
    oa_col: Optional[str] = None  # Other Allele
    eaf_col: Optional[str] = None  # Effect Allele Frequency

    def is_valid(self) -> bool:
        """检查是否包含最基本的必需列"""
        return all([self.snp_col, self.chr_col, self.bp_col, self.p_col])

    def to_dict(self) -> Dict:
        return {
            "SNP": self.snp_col,
            "CHR": self.chr_col,
            "BP": self.bp_col,
            "P": self.p_col,
            "BETA": self.beta_col,
            "SE": self.se_col,
            "EA": self.ea_col,
            "OA": self.oa_col,
            "EAF": self.eaf_col
        }

@dataclass
class DataPreview:
    """数据预览结果"""
    file_path: str
    file_format: FileFormat
    total_columns: int
    column_names: List[str]
    sample_rows: List[Dict]
    column_mapping: ColumnMapping
    file_size_mb: float
    estimated_rows: Optional[int]
    warnings: List[str]

    def to_dict(self) -> Dict:
        return {
            "file_path": self.file_path,
            "format": self.file_format.value,
            "columns": self.total_columns,
            "column_names": self.column_names,
            "sample_data": self.sample_rows,
            "detected_mapping": self.column_mapping.to_dict(),
            "file_size": f"{self.file_size_mb:.2f} MB",
            "estimated_rows": f"{self.estimated_rows:,}" if self.estimated_rows else "Unknown",
            "warnings": self.warnings
        }

class DataInspector:
    """GWAS 数据文件检查工具"""

    # 常见列名模式（正则表达式，不区分大小写）
    COLUMN_PATTERNS = {
        'snp': [r'^snp$', r'^rsid$', r'^rs$', r'^variant_id$', r'^marker$', r'^id$'],
        'chr': [r'^chr$', r'^chromosome$', r'^chrom$', r'^#chr$'],
        'bp': [r'^bp$', r'^pos$', r'^position$', r'^base_pair_location$'],
        'p': [r'^p$', r'^pval$', r'^p_value$', r'^pvalue$', r'^p\.value$', r'^p_bolt_lmm_inf$'],
        'beta': [r'^beta$', r'^b$', r'^effect$', r'^effect_size$', r'^beta_bolt_lmm_inf$'],
        'se': [r'^se$', r'^stderr$', r'^standard_error$', r'^se_bolt_lmm_inf$'],
        'ea': [r'^ea$', r'^effect_allele$', r'^a1$', r'^allele1$', r'^alt$'],
        'oa': [r'^oa$', r'^other_allele$', r'^a2$', r'^allele2$', r'^ref$'],
        'eaf': [r'^eaf$', r'^freq$', r'^maf$', r'^effect_allele_freq$', r'^a1freq$']
    }

    def __init__(self):
        self.warnings = []

    def peek_data(
        self,
        file_path: str,
        num_lines: int = 5
    ) -> DataPreview:
        """
        读取文件头部并分析结构

        Args:
            file_path: 文件路径
            num_lines: 读取的数据行数（不包括表头）

        Returns:
            DataPreview 对象
        """
        self.warnings = []
        file_path_obj = Path(file_path)

        if not file_path_obj.exists():
            raise FileNotFoundError(f"文件不存在: {file_path}")

        # 检测文件格式
        file_format = self._detect_format(file_path_obj)

        # 获取文件大小
        file_size_mb = file_path_obj.stat().st_size / (1024 * 1024)

        # 读取文件内容
        headers, rows = self._read_file_content(file_path_obj, file_format, num_lines)

        # 智能列名映射
        column_mapping = self._map_columns(headers)

        # 估算总行数
        estimated_rows = self._estimate_total_rows(file_path_obj, file_size_mb, len(rows))

        # 生成预览对象
        preview = DataPreview(
            file_path=file_path,
            file_format=file_format,
            total_columns=len(headers),
            column_names=headers,
            sample_rows=rows,
            column_mapping=column_mapping,
            file_size_mb=file_size_mb,
            estimated_rows=estimated_rows,
            warnings=self.warnings
        )

        return preview

    def _detect_format(self, file_path: Path) -> FileFormat:
        """检测文件格式"""
        suffix = file_path.suffix.lower()

        if suffix == '.gz':
            # 检查去掉 .gz 后的后缀
            stem_suffix = file_path.stem.split('.')[-1].lower()
            if stem_suffix == 'vcf':
                return FileFormat.VCF
            elif stem_suffix == 'csv':
                return FileFormat.CSV
            else:
                return FileFormat.TSV

        elif suffix == '.vcf':
            return FileFormat.VCF
        elif suffix == '.csv':
            return FileFormat.CSV
        elif suffix in ['.tsv', '.txt']:
            return FileFormat.TSV
        else:
            self.warnings.append(f"未知文件格式: {suffix}，尝试按 TSV 处理")
            return FileFormat.TSV

    def _read_file_content(
        self,
        file_path: Path,
        file_format: FileFormat,
        num_lines: int
    ) -> Tuple[List[str], List[Dict]]:
        """读取文件内容"""

        # 判断是否为压缩文件
        is_gzipped = file_path.suffix.lower() == '.gz'

        # 选择打开方式
        open_func = gzip.open if is_gzipped else open
        mode = 'rt' if is_gzipped else 'r'

        try:
            with open_func(file_path, mode, encoding='utf-8', errors='replace') as f:
                # 跳过 VCF 的元信息行
                if file_format == FileFormat.VCF:
                    for line in f:
                        if not line.startswith('##'):
                            header_line = line.strip()
                            break
                else:
                    header_line = f.readline().strip()

                # 解析表头
                delimiter = self._detect_delimiter(header_line, file_format)
                headers = header_line.lstrip('#').split(delimiter)
                headers = [h.strip() for h in headers]

                # 读取数据行
                rows = []
                for i, line in enumerate(f):
                    if i >= num_lines:
                        break

                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    values = line.split(delimiter)
                    row_dict = {headers[j]: values[j].strip() if j < len(values) else ''
                                for j in range(len(headers))}
                    rows.append(row_dict)

                return headers, rows

        except UnicodeDecodeError:
            self.warnings.append("文件编码异常，尝试使用 latin-1 编码")
            # 重试使用 latin-1 编码
            with open_func(file_path, mode, encoding='latin-1') as f:
                header_line = f.readline().strip()
                delimiter = self._detect_delimiter(header_line, file_format)
                headers = header_line.lstrip('#').split(delimiter)
                rows = []
                for i, line in enumerate(f):
                    if i >= num_lines:
                        break
                    values = line.strip().split(delimiter)
                    row_dict = {headers[j]: values[j] if j < len(values) else ''
                                for j in range(len(headers))}
                    rows.append(row_dict)
                return headers, rows

    def _detect_delimiter(self, header_line: str, file_format: FileFormat) -> str:
        """检测分隔符"""
        if file_format == FileFormat.CSV:
            return ','
        elif file_format == FileFormat.VCF:
            return '\t'
        else:
            # 自动检测 TSV 或空格分隔
            tab_count = header_line.count('\t')
            space_count = header_line.count(' ')

            if tab_count > space_count:
                return '\t'
            elif space_count > 0:
                self.warnings.append("检测到空格分隔符，可能导致解析问题")
                return ' '
            else:
                return '\t'

    def _map_columns(self, headers: List[str]) -> ColumnMapping:
        """智能映射列名到标准字段"""
        mapping = ColumnMapping()

        for header in headers:
            header_lower = header.lower().strip()

            # 尝试匹配每个标准字段
            for field, patterns in self.COLUMN_PATTERNS.items():
                for pattern in patterns:
                    if re.match(pattern, header_lower):
                        setattr(mapping, f"{field}_col", header)
                        break

        # 检查必需字段
        if not mapping.is_valid():
            missing = []
            if not mapping.snp_col:
                missing.append("SNP/rsID")
            if not mapping.chr_col:
                missing.append("Chromosome")
            if not mapping.bp_col:
                missing.append("Position")
            if not mapping.p_col:
                missing.append("P-value")

            self.warnings.append(f"⚠️ 缺少关键列: {', '.join(missing)}")

        return mapping

    def _estimate_total_rows(
        self,
        file_path: Path,
        file_size_mb: float,
        sample_rows: int
    ) -> Optional[int]:
        """估算文件总行数"""
        if sample_rows == 0:
            return None

        # 粗略估算：假设每行大小相似
        # 这只是一个粗略估计，实际可能有偏差
        avg_row_size_kb = (file_size_mb * 1024) / 1000  # 假设有1000行样本
        estimated = int((file_size_mb * 1024) / (avg_row_size_kb / sample_rows))

        return estimated

    def validate_gwas_format(self, file_path: str) -> Dict:
        """
        验证 GWAS 文件格式是否符合标准

        Returns:
            验证结果字典
        """
        try:
            preview = self.peek_data(file_path, num_lines=10)

            validation_result = {
                "is_valid": preview.column_mapping.is_valid(),
                "file_path": file_path,
                "detected_columns": preview.column_mapping.to_dict(),
                "missing_required": [],
                "warnings": preview.warnings,
                "recommendations": []
            }

            # 检查必需列
            if not preview.column_mapping.snp_col:
                validation_result["missing_required"].append("SNP identifier")
            if not preview.column_mapping.chr_col:
                validation_result["missing_required"].append("Chromosome")
            if not preview.column_mapping.bp_col:
                validation_result["missing_required"].append("Base pair position")
            if not preview.column_mapping.p_col:
                validation_result["missing_required"].append("P-value")

            # 生成建议
            if not preview.column_mapping.beta_col:
                validation_result["recommendations"].append(
                    "建议包含 BETA 列以进行效应量分析"
                )
            if not preview.column_mapping.se_col:
                validation_result["recommendations"].append(
                    "建议包含 SE 列以进行 meta 分析"
                )
            if not preview.column_mapping.ea_col or not preview.column_mapping.oa_col:
                validation_result["recommendations"].append(
                    "建议包含等位基因信息（EA/OA）以避免方向性错误"
                )

            return validation_result

        except Exception as e:
            return {
                "is_valid": False,
                "error": str(e),
                "file_path": file_path
            }


# Function Schemas for LLM
PEEK_DATA_SCHEMA = {
    "name": "peek_data_headers",
    "description": "读取 GWAS 数据文件的前几行，提取列名和样本数据，并智能识别标准字段（SNP, CHR, BP, P, BETA等）。在编写数据处理代码前必须先调用此工具。",
    "parameters": {
        "type": "object",
        "properties": {
            "file_path": {
                "type": "string",
                "description": "GWAS 数据文件的本地绝对路径，通常由 submit_download_task 返回。"
            },
            "num_lines": {
                "type": "integer",
                "description": "需要读取的数据行数（不包括表头），默认为 5 行。建议不超过 10 行以节省上下文。",
                "default": 5
            }
        },
        "required": ["file_path"]
    }
}

VALIDATE_FORMAT_SCHEMA = {
    "name": "validate_gwas_format",
    "description": "验证 GWAS 文件是否包含必需的标准列（SNP, CHR, BP, P），并提供格式改进建议。在开始分析前建议先调用此工具。",
    "parameters": {
        "type": "object",
        "properties": {
            "file_path": {
                "type": "string",
                "description": "需要验证的 GWAS 文件路径。"
            }
        },
        "required": ["file_path"]
    }
}


# 测试代码
if __name__ == "__main__":
    inspector = DataInspector()

    # 创建测试文件
    test_file = Path("./test_gwas.tsv")
    with open(test_file, 'w') as f:
        f.write("SNP\tCHR\tBP\tP\tBETA\tSE\tEA\tOA\n")
        f.write("rs123\t1\t1000\t0.001\t0.05\t0.01\tA\tG\n")
        f.write("rs456\t2\t2000\t0.05\t0.03\t0.02\tT\tC\n")
        f.write("rs789\t3\t3000\t1e-8\t0.1\t0.015\tG\tA\n")

    print("=== 测试数据预览 ===")
    preview = inspector.peek_data(str(test_file), num_lines=3)
    print(f"\n文件格式: {preview.file_formatvalue}")
    print(f"列数: {preview.total_columns}")
    print(f"列名: {preview.column_names}")
    print(f"\n检测到的标准列映射:")
    for key, value in preview.column_mapping.to_dict().items():
        print(f"  {key}: {value}")

    print(f"\n样本数据:")
    for i, row in enumerate(preview.sample_rows, 1):
        print(f"  行 {i}: {row}")

    print("\n=== 测试格式验证 ===")
    validation = inspector.validate_gwas_format(str(test_file))
    print(f"格式有效: {validation['is_valid']}")
    if validation['missing_required']:
        print(f"缺少必需列: {validation['missing_required']}")
    if validation['recommendations']:
        print("建议:")
        for rec in validation['recommendations']:
            print(f"  - {rec}")

    # 清理测试文件
    test_file.unlink()
