# code_executor.py
import sys
import io
import traceback
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from datetime import datetime
import json
import warnings
warnings.filterwarnings('ignore')

@dataclass
class ExecutionResult:
    """代码执行结果"""
    success: bool
    stdout: str
    stderr: str
    return_value: Any
    execution_time: float
    generated_files: List[str]
    error_traceback: Optional[str]
    variables_snapshot: Dict[str, str]

    def to_dict(self) -> Dict:
        return {
            "success": self.success,
            "stdout": self.stdout,
            "stderr": self.stderr,
            "return_value": str(self.return_value) if self.return_value else None,
            "execution_time": f"{self.execution_time:.3f}s",
            "generated_files": self.generated_files,
            "error": self.error_traceback,
            "variables": self.variables_snapshot
        }

class SafeCodeExecutor:
    """安全的代码执行器"""

    def __init__(self, output_dir: str = "./analysis_output"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 全局命名空间（用于保持变量状态）
        self.global_namespace = {
            'pd': pd,
            'np': np,
            'plt': plt,
            '__builtins__': __builtins__
        }
        
        # 执行历史
        self.execution_history = []

    def execute(
        self,
        code: str,
        description: str = "",
        timeout: int = 300
    ) -> ExecutionResult:
        """
        执行 Python 代码

        Args:
            code: 要执行的 Python 代码
            description: 代码描述（用于日志）
            timeout: 超时时间（秒）

        Returns:
            ExecutionResult 对象
        """
        start_time = datetime.now()
        generated_files = []
        
        # 捕获标准输出和错误
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()
        
        try:
            sys.stdout = stdout_capture
            sys.stderr = stderr_capture
            
            # 记录执行前的文件列表
            existing_files = set(self.output_dir.glob('*'))
            
            # 执行代码
            exec_result = None
            try:
                # 编译代码
                compiled_code = compile(code, '<string>', 'exec')
                
                # 执行
                exec(compiled_code, self.global_namespace)
                
                # 如果代码中有返回值表达式，尝试获取
                if 'result' in self.global_namespace:
                    exec_result = self.global_namespace['result']
                
                success = True
                error_traceback = None
                
            except Exception as e:
                success = False
                error_traceback = traceback.format_exc()
                print(f"执行错误: {str(e)}", file=sys.stderr)
            
            # 保存所有打开的图表
            if plt.get_fignums():
                for fig_num in plt.get_fignums():
                    fig = plt.figure(fig_num)
                    output_path = self.output_dir / f"plot_{fig_num}_{start_time.strftime('%Y%m%d_%H%M%S')}.png"
                    fig.savefig(output_path, dpi=300, bbox_inches='tight')
                    generated_files.append(str(output_path))
                    print(f"📊 图表已保存: {output_path}")
                plt.close('all')
            
            # 检测新生成的文件
            new_files = set(self.output_dir.glob('*')) - existing_files
            for f in new_files:
                if str(f) not in generated_files:
                    generated_files.append(str(f))
            
            # 计算执行时间
            execution_time = (datetime.now() - start_time).total_seconds()
            
            # 获取变量快照
            variables_snapshot = self._get_variables_snapshot()
            
            # 创建结果对象
            result = ExecutionResult(
                success=success,
                stdout=stdout_capture.getvalue(),
                stderr=stderr_capture.getvalue(),
                return_value=exec_result,
                execution_time=execution_time,
                generated_files=generated_files,
                error_traceback=error_traceback,
                variables_snapshot=variables_snapshot
            )
            
            # 记录历史
            self.execution_history.append({
                'timestamp': start_time.isoformat(),
                'description': description,
                'success': success,
                'execution_time': execution_time
            })
            
            return result
            
        finally:
            # 恢复标准输出
            sys.stdout = old_stdout
            sys.stderr = old_stderr

    def _get_variables_snapshot(self) -> Dict[str, str]:
        """获取当前命名空间中的关键变量"""
        snapshot = {}
        
        for name, value in self.global_namespace.items():
            # 跳过内置对象和模块
            if name.startswith('_') or name in ['pd', 'np', 'plt']:
                continue
            
            try:
                # 对于 DataFrame，显示形状和列名
                if isinstance(value, pd.DataFrame):
                    snapshot[name] = f"DataFrame(shape={value.shape}, columns={list(value.columns[:5])})"
                
                # 对于 Series
                elif isinstance(value, pd.Series):
                    snapshot[name] = f"Series(length={len(value)}, dtype={value.dtype})"
                
                # 对于数组
                elif isinstance(value, np.ndarray):
                    snapshot[name] = f"ndarray(shape={value.shape}, dtype={value.dtype})"
                
                # 对于列表和字典（限制长度）
                elif isinstance(value, (list, dict)):
                    snapshot[name] = f"{type(value).__name__}(length={len(value)})"
                
                # 其他简单类型
                elif isinstance(value, (int, float, str, bool)):
                    snapshot[name] = str(value)[:100]  # 限制字符串长度
                    
            except Exception:
                snapshot[name] = f"<{type(value).__name__}>"
        
        return snapshot

    def get_variable(self, var_name: str) -> Any:
        """获取命名空间中的变量"""
        return self.global_namespace.get(var_name)

    def clear_namespace(self):
        """清空命名空间（保留基础库）"""
        self.global_namespace = {
            'pd': pd,
            'np': np,
            'plt': plt,
            '__builtins__': __builtins__
        }
        print("✅ 命名空间已清空")


class GWASAnalysisHelper:
    """GWAS 分析辅助工具"""

    def __init__(self, executor: SafeCodeExecutor):
        self.executor = executor

    def load_gwas_data(
        self,
        file_path: str,
        column_mapping: Dict[str, str],
        nrows: Optional[int] = None
    ) -> ExecutionResult:
        """
        加载 GWAS 数据到 DataFrame

        Args:
            file_path: 文件路径
            column_mapping: 列名映射字典
            nrows: 读取行数限制

        Returns:
            ExecutionResult
        """
        code = f"""
import pandas as pd
import numpy as np

# 读取数据
print("🔄 正在加载数据...")
df = pd.read_csv(
    '{file_path}',
    sep='\\t',
    nrows={nrows if nrows else 'None'},
    low_memory=False
)

print(f"✅ 数据加载完成: {{df.shape[0]:,}} 行 × {{df.shape[1]}} 列")

# 重命名列
column_mapping = {column_mapping}
df = df.rename(columns=column_mapping)

# 数据类型转换
if 'CHR' in df.columns:
    df['CHR'] = df['CHR'].astype(str).str.replace('chr', '', case=False)
    df['CHR'] = pd.to_numeric(df['CHR'], errors='coerce')

if 'BP' in df.columns:
    df['BP'] = pd.to_numeric(df['BP'], errors='coerce')

if 'P' in df.columns:
    df['P'] = pd.to_numeric(df['P'], errors='coerce')

# 移除缺失值
initial_rows = len(df)
df = df.dropna(subset=['CHR', 'BP', 'P'])
removed_rows = initial_rows - len(df)

if removed_rows > 0:
    print(f"⚠️ 移除了 {{removed_rows:,}} 行缺失数据")

# 显示数据摘要
print("\\n数据摘要:")
print(df.head())
print(f"\\nP 值范围: {{df['P'].min():.2e}} - {{df['P'].max():.2e}}")
print(f"染色体: {{sorted(df['CHR'].unique())}}")

result = df
"""
        return self.executor.execute(code, description="加载 GWAS 数据")

    def manhattan_plot(
        self,
        significance_threshold: float = 5e-8,
        suggestive_threshold: float = 1e-5
    ) -> ExecutionResult:
        """
        生成曼哈顿图

        Args:
            significance_threshold: 显著性阈值
            suggestive_threshold: 提示性阈值

        Returns:
            ExecutionResult
        """
        code = f"""
import matplotlib.pyplot as plt
import numpy as np

# 检查数据是否存在
if 'df' not in globals():
    raise ValueError("请先使用 load_gwas_data 加载数据")

# 准备数据
plot_df = df[['CHR', 'BP', 'P']].copy()
plot_df = plot_df[plot_df['P'] > 0]  # 移除 P=0 的点
plot_df['-log10P'] = -np.log10(plot_df['P'])

# 计算每个染色体的 x 轴位置
plot_df['x_pos'] = 0
chr_centers = {{}}
x_offset = 0

for chrom in sorted(plot_df['CHR'].unique()):
    chr_data = plot_df[plot_df['CHR'] == chrom]
    chr_length = chr_data['BP'].max() - chr_data['BP'].min()
    
    plot_df.loc[plot_df['CHR'] == chrom, 'x_pos'] = chr_data['BP'] + x_offset
    chr_centers[chrom] = x_offset + chr_length / 2
    x_offset += chr_length + 1e7  # 染色体间隔

# 创建图表
fig, ax = plt.subplots(figsize=(16, 6))

# 按染色体着色
colors = ['#1f77b4', '#ff7f0e']
for i, chrom in enumerate(sorted(plot_df['CHR'].unique())):
    chr_data = plot_df[plot_df['CHR'] == chrom]
    ax.scatter(
        chr_data['x_pos'],
        chr_data['-log10P'],
        c=colors[int(chrom) % 2],
        s=5,
        alpha=0.7,
        label=f'Chr {{int(chrom)}}' if i < 2 else ''
    )

# 添加阈值线
sig_line = -np.log10({significance_threshold})
sug_line = -np.log10({suggestive_threshold})

ax.axhline(y=sig_line, color='red', linestyle='--', linewidth=1, label=f'Genome-wide significance (P={significance_threshold:.0e})')
ax.axhline(y=sug_line, color='blue', linestyle='--', linewidth=1, label=f'Suggestive (P={suggestive_threshold:.0e})')

# 设置 x 轴刻度
ax.set_xticks([chr_centers[c] for c in sorted(chr_centers.keys())])
ax.set_xticklabels([str(int(c)) for c in sorted(chr_centers.keys())])

# 标签和标题
ax.set_xlabel('Chromosome', fontsize=12, fontweight='bold')
ax.set_ylabel('-log10(P)', fontsize=12, fontweight='bold')
ax.set_title('Manhattan Plot', fontsize=14, fontweight='bold')
ax.legend(loc='upper right', fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()

# 统计显著位点
sig_snps = df[df['P'] < {significance_threshold}]
print(f"\\n📊 显著位点统计:")
print(f"  基因组显著 (P < {significance_threshold:.0e}): {{len(sig_snps):,}} 个")
print(f"  提示性显著 (P < {suggestive_threshold:.0e}): {{len(df[df['P'] < {suggestive_threshold}]):,}} 个")

if len(sig_snps) > 0:
    print(f"\\n  Top 5 显著位点:")
    top_snps = sig_snps.nsmallest(5, 'P')
    for idx, row in top_snps.iterrows():
        snp_id = row.get('SNP', 'N/A')
        print(f"    {{snp_id}}: Chr{{int(row['CHR'])}}:{{int(row['BP'])}} (P={{row['P']:.2e}})")
"""
        return self.executor.execute(code, description="生成曼哈顿图")

    def qq_plot(self) -> ExecutionResult:
        """生成 QQ 图"""
        code = """
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# 检查数据
if 'df' not in globals():
    raise ValueError("请先使用 load_gwas_data 加载数据")

# 准备数据
observed_p = df['P'].dropna().sort_values()
n = len(observed_p)

# 计算期望 P 值
expected_p = np.arange(1, n + 1) / (n + 1)

# 转换为 -log10
observed_log = -np.log10(observed_p)
expected_log = -np.log10(expected_p)

# 计算 lambda (基因组膨胀因子)
chi2_obs = stats.chi2.ppf(1 - observed_p, df=1)
lambda_gc = np.median(chi2_obs) / stats.chi2.ppf(0.5, df=1)

# 创建图表
fig, ax = plt.subplots(figsize=(8, 8))

# 绘制 QQ 图
ax.scatter(expected_log, observed_log, s=10, alpha=0.6, color='#1f77b4')

# 添加对角线
max_val = max(expected_log.max(), observed_log.max())
ax.plot([0, max_val], [0, max_val], 'r--', linewidth=2, label='Expected')

# 添加 lambda 信息
ax.text(
    0.05, 0.95,
    f'λ_GC = {lambda_gc:.3f}',
    transform=ax.transAxes,
    fontsize=12,
    verticalalignment='top',
    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5)
)

# 标签和标题
ax.set_xlabel('Expected -log10(P)', fontsize=12, fontweight='bold')
ax.set_ylabel('Observed -log10(P)', fontsize=12, fontweight='bold')
ax.set_title('QQ Plot', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()

print(f"\\n📊 QQ 图统计:")
print(f"  基因组膨胀因子 (λ_GC): {lambda_gc:.3f}")

if lambda_gc > 11:
    print(f"  ⚠️ λ_GC > 1.1，可能存在群体分层或其他混杂因素")
elif lambda_gc < 0.9:
    print(f"  ⚠️ λ_GC < 0.9，可能存在过度校正")
else:
    print(f"  ✅ λ_GC 在正常范围内")
"""
        return self.executor.execute(code, description="生成 QQ 图")

    def regional_plot(
        self,
        chromosome: int,
        start_bp: int,
        end_bp: int,
        lead_snp: Optional[str] = None
    ) -> ExecutionResult:
        """
        生成区域关联图

        Args:
            chromosome: 染色体编号
            start_bp: 起始位置
            end_bp: 结束位置
            lead_snp: 先导 SNP（可选）

        Returns:
            ExecutionResult
        """
        code = f"""
import matplotlib.pyplot as plt
import numpy as np

# 检查数据
if 'df' not in globals():
    raise ValueError("请先使用 load_gwas_data 加载数据")

# 提取区域数据
region_df = df[
    (df['CHR'] == {chromosome}) &
    (df['BP'] >= {start_bp}) &
    (df['BP'] <= {end_bp})
].copy()

if len(region_df) == 0:
    raise ValueError(f"在 Chr{chromosome}:{start_bp}-{end_bp} 未找到数据")

region_df['-log10P'] = -np.log10(region_df['P'])

# 创建图表
fig, ax = plt.subplots(figsize=(12, 6))

# 绘制散点
scatter = ax.scatter(
    region_df['BP'] / 1e6,  # 转换为 Mb
    region_df['-log10P'],
    c=region_df['-log10P'],
    cmap='YlOrRd',
    s=30,
    alpha=0.7,
    edgecolors='black',
    linewidth=0.5
)

# 添加颜色条
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('-log10(P)', fontsize=10)

# 标记先导 SNP
lead_snp_id = '{lead_snp}' if '{lead_snp}' != 'None' else None
if lead_snp_id and 'SNP' in region_df.columns:
    lead_data = region_df[region_df['SNP'] == lead_snp_id]
    if not lead_data.empty:
        ax.scatter(
            lead_data['BP'] / 1e6,
            lead_data['-log10P'],
            s=200,
            marker='D',
            c='purple',
            edgecolors='black',
            linewidth=2,
            label=f'Lead SNP: {{lead_snp_id}}',
            zorder=5
        )

# 添加显著性阈值线
ax.axhline(y=-np.log10(5e-8), color='red', linestyle='--', linewidth=1, label='P=5e-8')

# 标签和标题
ax.set_xlabel(f'Position on Chromosome {chromosome} (Mb)', fontsize=12, fontweight='bold')
ax.set_ylabel('-log10(P)', fontsize=12, fontweight='bold')
ax.set_title(f'Regional Association Plot: Chr{chromosome}:{start_bp/1e6:.2f}-{end_bp/1e6:.2f} Mb', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()

print(f"\\n📊 区域统计:")
print(f"  区域: Chr{chromosome}:{start_bp:,}-{end_bp:,}")
print(f"  SNP 数量: {{len(region_df):,}}")
print(f"  最小 P 值: {{region_df['P'].min():.2e}}")

# 显示 Top 5 SNPs
top_snps = region_df.nsmallest(5, 'P')
print(f"\\n  Top 5 SNPs:")
for idx, row in top_snps.iterrows():
    snp_id = row.get('SNP', 'N/A')
    print(f"    {{snp_id}}: {{int(row['BP']):,}} bp (P={{row['P']:.2e}})")
"""
        return self.executor.execute(code, description="生成区域关联图")


# Function Schema for LLM
EXECUTE_CODE_SCHEMA = {
    "name": "execute_analysis_code",
    "description": "在隔离的 Python 环境中执行 GWAS 分析代码。支持 pandas、numpy、matplotlib 等库。代码执行结果会被捕获并返回，生成的图表会自动保存。",
    "parameters": {
        "type": "object",
        "properties": {
            "code": {
                "type": "string",
                "description": "要执行的 Python 代码。可以使用 df（DataFrame）、pd（pandas）、np（numpy）、plt（matplotlib）等预定义变量。"
            },
            "description": {
                "type": "string",
                "description": "代码功能描述，用于日志记录。",
                "default": ""
            }
        },
        "required": ["code"]
    }
}


# 测试代码
if __name__ == "__main__":
    # 创建执行器
    executor = SafeCodeExecutor(output_dir="./test_output")
    helper = GWASAnalysisHelper(executor)

    # 创建测试数据
    test_data_path = Path("./test_gwas_large.tsv")
    
    print("=== 创建测试数据 ===")
    np.random.seed(42)
    n_snps = 10000
    
    test_df = pd.DataFrame({
        'SNP': [f'rs{i}' for i in range(n_snps)],
        'CHR': np.random.choice(range(1, 23), n_snps),
        'BP': np.random.randint(1000000, 100000000, n_snps),
        'P': np.random.beta(0.1, 10, n_snps),  # 生成偏向小值的 P 值
        'BETA': np.random.normal(0, 0.1, n_snps),
        'SE': np.random.uniform(0.01, 0.05, n_snps)
    })
    
    # 添加一些显著信号
    test_df.loc[:50, 'P'] = np.random.uniform(1e-10, 1e-7, 51)
    
    test_df.to_csv(test_data_path, sep='\t', index=False)
    print(f"✅ 测试数据已创建: {test_data_path}")

    # 测试 1: 加载数据
    print("\n=== 测试 1: 加载数据 ===")
    result = helper.load_gwas_data(
        str(test_data_path),
        column_mapping={'SNP': 'SNP', 'CHR': 'CHR', 'BP': 'BP', 'P': 'P'}
    )
    print(f"执行状态: {'成功' if result.success else '失败'}")
    print(f"执行时间: {result.execution_time:.3f}s")
    print(f"输出:\n{result.stdout}")

    # 测试 2: 曼哈顿图
    print("\n=== 测试 2: 生成曼哈顿图 ===")
    result = helper.manhattan_plot()
    print(f"执行状态: {'成功' if result.success else '失败'}")
    print(f"生成文件: {result.generated_files}")
    print(f"输出:\n{result.stdout}")

    # 测试 3: QQ 图
    print("\n=== 测试 3: 生成 QQ 图 ===")
    result = helper.qq_plot()
    print(f"执行状态: {'成功' if result.success else '失败'}")
    print(f"生成文件: {result.generated_files}")
    print(f"输出:\n{result.stdout}")

    # 测试 4: 自定义代码
    print("\n=== 测试 4: 执行自定义分析 ===")
    custom_code = """
# 计算每个染色体的显著 SNP 数量
sig_threshold = 5e-8
sig_by_chr = df[df['P'] < sig_threshold].groupby('CHR').size()

print("\\n各染色体显著 SNP 数量:")
for chr_num, count in sig_by_chr.items():
    print(f"  Chr {int(chr_num)}: {count} 个")

result = sig_by_chr
"""
    result = executor.execute(custom_code, description="统计显著 SNP")
    print(f"执行状态: {'成功' if result.success else '失败'}")
    print(f"输出:\n{result.stdout}")

    # 清理测试文件
    test_data_path.unlink()
    print("\n✅ 测试完成")
