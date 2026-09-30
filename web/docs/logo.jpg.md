# `static/logo.jpg`：品牌图片

页面左上角 DPTechnology／深势科技图片，由 `index.html` 引用，尺寸由 `styles.css` 的 `.brand img` 控制。服务器部署时需将图片随静态资源一起发布，并确认有权在目标环境使用该品牌图片。若改用透明图或 SVG，更新 HTML 路径、缓存版本和替代文本，检查深浅背景下的可读性。
