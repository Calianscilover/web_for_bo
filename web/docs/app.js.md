# `static/app.js`：前端流程控制

`state` 保存当前 `designId`、`runId`、目标模式、候选页码和轮次。`call`/`post` 是统一 API 请求入口；`configFromForm` 把配方表单转换为步骤 1 的 JSON；`createDesign`、`loadCandidates`、`uploadExperiment` 处理候选和初始实验；`startRun`、`pollRun`、`loadRecommendation` 负责异步运行及产物预览；`submitFeedback`、`simulate`、`nextRound` 处理连续轮次。

前端可编辑推荐性能值，提交时按 `sample_id` 和目标列名构造 `measurements`；也能上传填好的推荐 CSV。回传成功后刷新推荐和可视化。`displayTable` 只显示部分便于阅读的列，完整列仍在下载 CSV 中。表格内容通过 `escapeHtml` 处理后插入页面。

迁服务器时，首先把 `api='/api/v1'` 变成构建环境配置；其次考虑任务断线重连、身份认证与权限、长耗时任务状态展示。不要只在浏览器内保存 run ID；当前后端目录已能列出设计与运行，后期应使用用户/项目关联。
