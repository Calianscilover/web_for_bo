# `static/visualization.css`：回传和可视化样式

负责推荐表下方的实验回传表单、按钮禁用态、两栏图表、Pareto 详情、预测列表、历史实验表和来源标记。关键类包括 `.feedback-entry`、`.viz-grid`、`.chart`、`.pareto-detail`、`.prediction-list`、`.source-badge`。

图表使用 `visualization.js` 生成 SVG，本文件控制容器、轴线、网格和自适应尺寸。修改图表颜色时同步 `visualization.js` 的 `chartColors`，以保持初始/真实/模拟来源含义一致。
