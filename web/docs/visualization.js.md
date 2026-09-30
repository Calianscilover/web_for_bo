# `static/visualization.js`：图表渲染

入口 `window.renderVisualization(data)` 接收 `GET /api/v1/optimization-runs/{run_id}/visualization` 返回的 JSON。`renderFitChart` 绘制训练样本「实测 vs 后验均值」散点；`renderProgressChart` 绘制单目标累计最佳实测性能或双目标 Pareto 超体积；`renderParetoChart` 绘制当前已测 Pareto 前沿；`renderPredictions` 和 `renderHistory` 展示本轮预测与历史实验。模拟、初始、真实回传使用不同颜色。

数据契约包括 `mode`、`targets`、`directions`、`model_round`、`fit`、`progress`、`pareto`、`measurements`、`recommendations`、`counts`，双目标含 `reference_point`。JSON 由 `visualize_bo.py` 计算，前端只绘图。训练拟合图不能当作留出集泛化精度；双目标超体积需要固定参考点和明确优化方向。

未来若改用图表库，保持 API 数据契约并替换本文件即可。若新增不确定度或跨轮模型比较，应先在步骤 4 和 `visualize_bo.py` 定义计算，再加前端表现。
