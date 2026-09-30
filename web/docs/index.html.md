# `static/index.html`：页面结构

单页工作台按「配方空间 → 候选池 → 实验数据 → 优化设置 → 推荐配方 → 可视化」排列。顶部使用 `logo.jpg`，页面通过 `styles.css`、`visualization.css`、`visualization.js`、`app.js` 加载资源。

关键 DOM 区域为 `#space`、`#pool`、`#experiment`、`#optimization`、`#recommendations`、`#visualization`。`app.js` 直接通过元素 ID 绑定事件，修改 ID 必须同步 JS。设计表单构成步骤 1 的 JSON；实验模板/上传区对应步骤 3；推荐、回传和图表区对应步骤 4 的输出。

迁服务器后页面可以仍由 Flask 静态文件服务，也可以由 Nginx/CDN 独立提供；独立部署时须配置 `/api/v1` 请求目标和同源策略，并保持下载链接指向有权限的 API。
