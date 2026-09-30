# `static/visualization.js`：图表渲染

入口 `window.renderVisualization(data)` 接收 `GET /api/v1/optimization-runs/{run_id}/visualization` 返回的 JSON。`renderFitChart` 绘制训练样本「实测 vs 后验均值」散点；`renderProgressChart` 绘制单目标累计最佳实测性能或多目标 Pareto 超体积；`renderParetoChart` 绘制全部已测配方（淡色点）与当前 Pareto 前沿（实心点，两个目标时以折线相连；3–4 个目标时由 `paretoAxes` 提供横纵轴选择，图为按全部目标判定的前沿在两个目标上的投影），最新回传轮次的配方加虚线圈，图下 `#pareto-note` 说明该轮新增多少条、有几条进入前沿；`renderPredictions` 和 `renderHistory` 展示本轮预测与历史实验。模拟、初始、真实回传使用不同颜色。`chartRange`/`chartTicks` 把坐标范围和刻度对齐到 1、2、2.5、5 × 10ⁿ 的整齐步长；进展图只有一个点时显示回传提示。`compositionHtml` 渲染配方组成标签，由历史表和 `app.js` 的候选池/推荐表共用。

数据契约包括 `mode`、`targets`、`directions`、`model_round`、`fit`、`progress`、`pareto`、`measurements`、`recommendations`、`counts`，多目标含 `reference_point`。JSON 由 `visualize_bo.py` 计算，前端只绘图。训练拟合图不能当作留出集泛化精度；多目标超体积需要固定参考点和明确优化方向。

未来若改用图表库，保持 API 数据契约并替换本文件即可。若新增不确定度或跨轮模型比较，应先在步骤 4 和 `visualize_bo.py` 定义计算，再加前端表现。
