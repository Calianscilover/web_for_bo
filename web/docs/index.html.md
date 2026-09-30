# `static/index.html`：页面结构

单页工作台按「配方空间 → 候选池 → 实验数据与优化目标 → 模型优化 → 推荐配方 → 可视化」排列。`.shell` 分为左侧固定流程栏 `.side-nav` 与内容区 `main`；流程栏 `#progress` 的六个链接（`data-step` 1–6）与卡片编号一一对应，既是唯一的页内导航，也显示各步完成状态。顶栏只保留 logo 与运行状态。页面通过 `styles.css`、`visualization.css`、`visualization.js`、`app.js` 加载资源。

关键 DOM 区域为 `#space`、`#pool`、`#experiment`、`#optimization`、`#recommendations`、`#visualization`。`app.js` 直接通过元素 ID 绑定事件，修改 ID 必须同步 JS。锂盐、溶剂与添加剂（`#additives`，可选，由 `#add-additive` 添加）用紧凑的 `.component-row` 行编辑，构成步骤 1 的 JSON。单/双目标切换、目标列名和方向集中在 `#experiment`（目标 2 行 `#target-2-row` 仅双目标显示），因为实验模板与上传都依赖目标列；实验上传只有一个入口 `#experiment-file`：`#template`（候选池模板，列出全部候选）与 `#recipe-template`（自有配方模板，含 EXAMPLE 示例）列完全相同，都随当前设计的组分变化；`#mapping-link` 下载对齐结果。`#optimization` 只显示只读摘要 `#target-summary`；「模型超参数设置」`#hyperparameters` 按高斯过程核、噪声与拟合、采集与候选评分、超体积参考点（仅双目标 `#hp-ref-group`）分组，噪声与参考点输入按目标名标注。推荐区的回传操作在 `#feedback-actions`，`#next-round-bar` 说明能否生成下一轮。

迁服务器后页面可以仍由 Flask 静态文件服务，也可以由 Nginx/CDN 独立提供；独立部署时须配置 `/api/v1` 请求目标和同源策略，并保持下载链接指向有权限的 API。
