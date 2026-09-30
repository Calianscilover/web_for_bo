# `static/app.js`：前端流程控制

`state` 保存当前 `designId`、`runId`、目标模式、候选页码、轮次，以及 `designReady`、`hasExperiment`、`feedbackComplete` 等进度标记；`updateProgress` 据此刷新左侧流程栏与卡片编号的完成/当前状态；`watchSections` 用 `IntersectionObserver` 给当前可见卡片对应的流程项加 `viewing`，点击流程项即滚动到对应卡片。`setMode` 切换单/双目标并刷新 `updateTargetSummary`；`uploadExperiment` 把两种模板填写后的 CSV 都提交到配方导入接口并刷新候选池，`afterExperimentImport` 在目标值完整的实验不少于 2 条时跳转到「模型优化」；`setExperimentLocked` 在设计已有运行时禁用上传；`additiveRow` 管理可选添加剂行；`updateNextRound` 根据运行的 `feedback` 状态说明为什么可以或不能生成下一轮。`call`/`post` 是统一 API 请求入口；`configFromForm` 把配方表单转换为步骤 1 的 JSON；`createDesign`、`loadCandidates`、`uploadExperiment` 处理候选和初始实验；`startRun`、`pollRun`、`loadRecommendation` 负责异步运行及产物预览；`submitFeedback`、`simulate`、`nextRound` 处理连续轮次。

前端可编辑推荐性能值，提交时按 `sample_id` 和目标列名构造 `measurements`；也能上传填好的推荐 CSV。回传成功后刷新推荐和可视化。`recipeTable` 把 `compound_i`/`mass_ratio_i` 合并成“组分 + 质量百分数”的组成列（只列非零组分，调用 `visualization.js` 的 `compositionHtml`），推荐表额外显示目标列，完整列仍在下载 CSV 中。表格内容通过 `escapeHtml` 处理后插入页面。

迁服务器时，首先把 `api='/api/v1'` 变成构建环境配置；其次考虑任务断线重连、身份认证与权限、长耗时任务状态展示。不要只在浏览器内保存 run ID；当前后端目录已能列出设计与运行，后期应使用用户/项目关联。
