# `web/test_integration.py`：流程验证

Flask `test_client` 与临时输出目录覆盖新设计撒点、候选转换、实验 CSV 上传、单目标/双目标/三目标首轮推荐（含目标数量与参考点长度校验）、真实/模拟实验回传、第二轮及第三轮结果、可视化和多目标超体积。运行：

```bash
python web/test_integration.py
```

在仓库根目录、已激活 `botorch` 环境（见 `environment.yml`）下运行。测试需要其中的 Flask、RDKit、PyTorch、BoTorch 等依赖，并会实际启动模型拟合子进程。它不访问或覆盖 `bo_test` 真实任务数据。迁服务器后可继续把此测试作为本地回归检查，并增加独立 worker、持久存储和重启恢复的部署级测试。
