# query_open_gwas_api.py
import requests
from typing import List, Dict, Optional
from dataclasses import dataclass
import json

@dataclass
class GWASDataset:
    """GWAS 数据集结构"""
    id: str
    trait: str
    sample_size: int
    ncase: Optional[int]
    ncontrol: Optional[int]
    population: str
    year: Optional[int]
    author: str
    consortium: Optional[str]
    pmid: Optional[str]
    
    def to_dict(self) -> Dict:
        return {
            "Dataset ID": self.id,
            "Trait": self.trait,
            "Sample Size": f"{self.sample_size:,}" if self.sample_size else "N/A",
            "Case/Control": f"{self.ncase}/{self.ncontrol}" if self.ncase and self.ncontrol else "Continuous trait",
            "Population": self.population,
            "Year": self.year or "N/A",
            "Author": self.author,
            "PMID": self.pmid or "N/A"
        }
class OpenGWASQueryTool:
    """OpenGWAS API 查询工具"""
    
    BASE_URL = "https://gwas.mrcieu.ac.uk/api"
    
    # 常见中文表型到英文的映射
    PHENOTYPE_TRANSLATION = {
        "2型糖尿病": "Type 2 Diabetes",
        "身高": "Height",
        "体重指数": "Body Mass Index",
        "精神分裂症": "Schizophrenia",
        "冠心病": "Coronary Artery Disease",
        "高血压": "Hypertension",
        "阿尔茨海默病": "Alzheimer's Disease",
        "帕金森病": "Parkinson's Disease",
        "哮喘": "Asthma",
        "类风湿关节炎": "Rheumatoid Arthritis",
        "海马体":"hippocampus"
    }
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'OpenClaw-GWAS-Agent/1.0'
        })
    
    def translate_phenotype(self, phenotype: str) -> str:
        """将中文表型翻译为英文"""
        return self.PHENOTYPE_TRANSLATION.get(phenotype, phenotype)
    
    def query(
        self, 
        phenotype: str, 
        population: str = "Any",
        top_k: int = 5
    ) -> List[GWASDataset]:
        """
        查询 OpenGWAS 数据库
        
        Args:
            phenotype: 表型名称
            population: 族裔筛选
            top_k: 返回结果数量
            
        Returns:
            GWASDataset 对象列表
        """
        # 翻译中文表型
        english_phenotype = self.translate_phenotype(phenotype)
        
        try:
            # 调用 OpenGWAS API
            response = self.session.get(
                f"{self.BASE_URL}/gwasinfo",
                params={"trait": english_phenotype},
                timeout=30
            )
            response.raise_for_status()
            
            data = response.json()
            
            if not data:
                print(f"⚠️ 未找到与 '{phenotype}' 相关的 GWAS 数据集")
                return []
            
            # 解析并过滤结果
            datasets = self._parse_results(data, population)
            
            # 排序：优先样本量大、近期发布
            datasets.sort(
                key=lambda x: (
                    x.sample_size if x.sample_size else 0,
                    x.year if x.year else 0
                ),
                reverse=True
            )
            
            return datasets[:top_k]
            
        except requests.exceptions.RequestException as e:
            print(f"API 请求失败: {str(e)}")
            return []
        except json.JSONDecodeError:
            print("API 返回数据格式错误")
            return []
    
    def _parse_results(self, data: List[Dict], population_filter: str) -> List[GWASDataset]:
        """解析 API 返回的原始数据"""
        datasets = []
        
        for item in data:
            # 族裔筛选
            pop = item.get("population", "Unknown")
            if population_filter != "Any" and population_filter.lower() not in pop.lower():
                continue
            
            dataset = GWASDataset(
                id=item.get("id", ""),
                trait=item.get("trait", ""),
                sample_size=item.get("sample_size", 0),
                ncase=item.get("ncase"),
                ncontrol=item.get("ncontrol"),
                population=pop,
                year=item.get("year"),
                author=item.get("author", "Unknown"),
                consortium=item.get("consortium"),
                pmid=item.get("pmid")
            )
            datasets.append(dataset)
        
        return datasets
    
    def format_results(self, datasets: List[GWASDataset]) -> str:
        """格式化输出结果"""
        if not datasets:
            return "未找到相关数据集"
        
        output = f"\n找到 {len(datasets)} 个相关的 GWAS 数据集：\n\n"
        
        for i, ds in enumerate(datasets, 1):
            output += f"【{i}】 {ds.trait}\n"
            output += f"  Dataset ID: {ds.id}\n"
            output += f"  样本量: {ds.sample_size:,}\n"
            
            if ds.ncase and ds.ncontrol:
                output += f"  病例/对照: {ds.ncase:,} / {ds.ncontrol:,}\n"
            
            output += f"  族裔: {ds.population}\n"
            output += f"  作者: {ds.author}"
            
            if ds.year:
                output += f" ({ds.year})"
            
            if ds.pmid:
                output += f"\n  PMID: {ds.pmid}"
            
            output += "\n\n"
        
        return output


# Function Schema for LLM
FUNCTION_SCHEMA = {
    "name": "query_open_gwas_api",
    "description": "通过 OpenGWAS 数据库检索与特定表型（Phenotype）或疾病相关的全基因组关联分析（GWAS）汇总数据集。支持中英文表型输入。返回数据集的ID、样本量、所属文献等信息。",
    "parameters": {
        "type": "object",
        "properties": {
            "phenotype": {
                "type": "string",
                "description": "表型名称或疾病名称，支持中文（如'2型糖尿病'）或英文（如'Type 2 Diabetes'）。"
            },
            "population": {
                "type": "string",
                "enum": ["European", "East Asian", "African", "Mixed", "Any"],
                "description": "筛选特定族裔的队列数据，默认填 'Any'。"
            },
            "top_k": {
                "type": "integer",
                "description": "返回最相关的结果数量，建议设置为 5。",
                "default": 5
            }
        },
        "required": ["phenotype"]
    }
}


# 测试代码
if __name__ == "__main__":
    tool = OpenGWASQueryTool()
    
    # 测试案例 1: 英文表型
    print("=== 测试 1: Type 2 Diabetes ===")
    results = tool.query("Type 2 Diabetes", population="European", top_k=3)
    print(tool.format_results(results))
    
    # 测试案例 2: 中文表型
    print("\n=== 测试 2: 2型糖尿病 ===")
    results = tool.query("2型糖尿病", top_k=3)
    print(tool.format_results(results))
    
    # 测试案例 3: 连续性状
    print("\n=== 测试 3: Body Mass Index ===")
    results = tool.query("Body Mass Index", population="East Asian", top_k=3)
    print(tool.format_results(results))
