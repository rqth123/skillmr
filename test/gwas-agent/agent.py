# gwas_agent/agent.py
"""
GWAS Analysis Agent - 主控制器
"""

from typing import List, Dict, Optional
from openclaw import Agent, Skill
from .skills import SKILL_DEFINITIONS


class GWASAgent(Agent):
    """
    GWAS 分析专用 Agent
    
    功能：
    1. 搜索和发现 GWAS 数据
    2. 下载和管理数据文件
    3. 数据质量检查和格式验证
    4. 执行统计分析和可视化
    """
    
    def __init__(
        self,
        name: str = "GWAS-Analyst",
        model: str = "claude-3-5-sonnet-20241022",
        **kwargs
    ):
        """
        初始化 GWAS Agent
        
        Args:
            name: Agent 名称
            model: 使用的 LLM 模型
            **kwargs: 其他 OpenClaw Agent 参数
        """
        
        # 构建系统提示词
        system_prompt = self._build_system_prompt()
        
        # 初始化父类
        super().__init__(
            name=name,
            model=model,
            system_prompt=system_prompt,
            **kwargs
        )
        
        # 注册技能
        self._register_skills()
    
    def _build_system_prompt(self) -> str:
        """构建系统提示词"""
        return """你是一个专业的 GWAS（全基因组关联分析）数据分析助手。

你的核心能力：
1. 搜索和发现 GWAS 数据集
2. 下载和管理大型基因组数据文件
3. 数据质量检查和格式验证
4. 执行统计分析（曼哈顿图、QQ图、区域关联图等）
5. 解释分析结果并提供生物学见解

工作流程：
1. 理解用户的研究问题
2. 搜索相关的 GWAS 数据
3. 下载并验证数据质量
4. 执行适当的统计分析
5. 生成可视化结果
6. 解释发现并提供建议

注意事项：
- 始终先验证数据格式再进行分析
- 对于大文件，使用异步下载
- 生成的图表会自动保存到本地
- 解释结果时要考虑多重检验校正
- 提醒用户注意群体分层等混杂因素
-不能有任何的编造和胡诌,一切以提供以及下载的数据为基础
你的回答应该：
- 专业但易懂
- 提供具体的代码和命令
- 解释统计概念
- 给出可操作的建议
"""
    
    def _register_skills(self):
        """注册所有技能"""
        for skill_def in SKILL_DEFINITIONS:
            skill = Skill(
                name=skill_def["name"],
                function=skill_def["function"],
                description=skill_def["description"],
                parameters=skill_def["parameters"]
            )
            self.add_skill(skill)
    
    def analyze_trait(
        self,
        trait_name: str,
        auto_download: bool = True,
        generate_plots: bool = True
    ) -> Dict:
        """
        一键分析特定表型的 GWAS 数据
        
        Args:
            trait_name: 表型名称
            auto_download: 是否自动下载数据
            generate_plots: 是否生成可视化
        
        Returns:
            分析结果字典
        """
        workflow = f"""
请帮我分析 '{trait_name}' 的 GWAS 数据：

1. 搜索相关的 GWAS 研究
2. {'下载最相关的数据集' if auto_download else '列出可用的数据集'}
3. 验证数据格式
4. {'生成曼哈顿图和 QQ 图' if generate_plots else '进行基础统计分析'}
5. 总结主要发现

请逐步执行并报告进度。
"""
        return self.run(workflow)


# 便捷函数
def create_gwas_agent(**kwargs) -> GWASAgent:
    """创建 GWAS Agent 实例"""
    return GWASAgent(**kwargs)
