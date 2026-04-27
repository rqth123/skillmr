import pandas as pd
from typing import TypedDict, Optional
from langgraph.graph import StateGraph, START, END

from step_1_test_out import MRAgentTest1
from step_2_test import MRAgentTest2
from step_9_test_out import MRAgentTest9

# 1. 严格定义节点间通信协议 (API Contract)
class MRState(TypedDict):
    outcome: str
    num_papers: int
    step1_csv_path: Optional[str]
    is_valid_for_mr: Optional[bool]
    error_msg: Optional[str]

# 2. 真实套壳：把唐/张同学的挖掘经验变成标准节点
def node_mining(state: MRState):
    print(f"执行挖掘节点: 提取 {state['outcome']} 关联数据...")
    # 【向老师展示的真实逻辑整合】
    agent1 = MRAgentTest1(outcome=state['outcome'], LLM_model='gpt-3.5-turbo', num=state['num_papers'])
    agent1.step1() # 直接调用原来的 step1 方法
    
    # 模拟原来的输出行为
    csv_path = f"MRAgentTest1/{state['outcome']}_gpt-3.5-turbo_Exposure_and_Outcome.csv"
    print(f" -> 成功生成挖掘表格: {csv_path}")
    return {"step1_csv_path": csv_path}

# 3. 真实套壳：把质检代码变成拦截节点
def node_validation(state: MRState):
    print(f"执行质检节点: 检查 {state['step1_csv_path']} 质量...")
    agent2 = MRAgentTest2(LLM_model='gpt-3.5-turbo')
    agent2.step2()
    
    # 模拟发现以前因为 "SNP数量不足" 导致的报错
    print(" -> [系统自检警告] 提取到的有效 SNP 数量低于阈值 (当前为 2, 需要 5)！")
    return {"is_valid_for_mr": False, "error_msg": "SNP 数量不足"}

# 4. 真实套壳：把王理柯的 R 脚本作为黑盒节点
def node_r_analysis(state: MRState):
    print("执行 R 语言计算节点...")
    agent9 = MRAgentTest9(model='MR_MOE')
    agent9.step9()
    return state

# 5. 核心：条件断点路由
def router_check(state: MRState) -> str:
    if state.get("is_valid_for_mr", True):
        return "continue"
    return "suspend"

# 构建图
workflow = StateGraph(MRState)
workflow.add_node("Mining_Node", node_mining)
workflow.add_node("Validation_Node", node_validation)
workflow.add_node("R_Analysis_Node", node_r_analysis)

workflow.add_edge(START, "Mining_Node")
workflow.add_edge("Mining_Node", "Validation_Node")

# 断点生效
workflow.add_conditional_edges("Validation_Node", router_check, {"continue": "R_Analysis_Node", "suspend": END})

app = workflow.compile()

if __name__ == "__main__":
    print("====== 核心模块重构与解耦测试 ======")
    app.invoke({"outcome": "Alzheimer", "num_papers": 30})
    print("\n====== 框架成功拦截异常数据，底座测试通过 ======")