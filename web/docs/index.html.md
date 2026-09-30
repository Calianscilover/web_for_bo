# `static/index.html`：页面结构

单页工作台按「配方空间 → 候选池 → 实验数据与优化目标 → 模型优化 → 推荐配方 → 可视化」排列，六张卡片的编号与标题右侧的进度条 `#progress` 一一对应（`data-step` 1–6）。顶部使用 `logo.jpg`，页面通过 `styles.css`、`visualization.css`、`visualization.js`、`app.js` 加载资源。

关键 DOM 区域为 `#space`、`#pool`、`#experiment`、`#optimization`、`#recommendations`、`#visualization`。`app.js` 直接通过元素 ID 绑定事件，修改 ID 必须同步 JS。锂盐与溶剂用紧凑的 `.component-row` 行编辑，构成步骤 1 的 JSON。单/双目标切换、目标列名和方向集中在 `#experiment`（目标 2 行 `#target-2-row` 仅双目标显示），因为实验模板与上传都依赖目标列；`#optimization` 只显示只读摘要 `#target-summary`。推荐区的回传操作在 `#feedback-actions`，`#next-round-bar` 说明能否生成下一轮。

迁服务器后页面可以仍由 Flask 静态文件服务，也可以由 Nginx/CDN 独立提供；独立部署时须配置 `/api/v1` 请求目标和同源策略，并保持下载链接指向有权限的 API。
