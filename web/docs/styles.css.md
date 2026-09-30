# `static/styles.css`：全局视觉样式

负责顶栏、logo 尺寸、页面网格、卡片、表单、候选/推荐表格、按钮、状态提示和响应式布局。`@media(max-width:1000px)` 将双列页面改为单列；`@media(max-width:600px)` 缩小移动端边距和字体。

后续品牌视觉调整以根部颜色和 `.brand`、`.topbar`、`.primary` 为入口。新增交互状态时优先沿用现有 `.notice`、`.disabled`、`.secondary` 风格，并检查键盘焦点、对比度和窄屏横向滚动。图表及实验回传专有样式放在 `visualization.css`。
