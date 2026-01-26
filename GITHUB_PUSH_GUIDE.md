# GitHub 推送指南

## 仓库信息

| 项目 | 值 |
|------|-----|
| **仓库地址** | https://github.com/cruan1991/SCONE-DNA |
| **用户名** | cruan1991 |
| **本地路径** | /Users/mac/Documents/SCONE-DNA |
| **认证方式** | GitHub CLI (gh) ✅ 已配置 |

---

## 日常更新流程（一键推送）

```bash
cd /Users/mac/Documents/SCONE-DNA
git add .
git commit -m "你的更新说明"
git push
```

**不需要输入密码或 token！** GitHub CLI 已配置自动认证。

---

## 详细步骤

### 1. 进入项目目录
```bash
cd /Users/mac/Documents/SCONE-DNA
```

### 2. 查看修改状态
```bash
git status
```

### 3. 添加修改的文件
```bash
# 添加所有修改
git add .

# 或只添加特定文件
git add 文件名.py
```

### 4. 提交修改
```bash
git commit -m "简短描述你的修改"
```

### 5. 推送到 GitHub
```bash
git push
```

---

## 常用命令速查

| 命令 | 说明 |
|------|------|
| `git status` | 查看当前状态 |
| `git add .` | 添加所有修改 |
| `git commit -m "msg"` | 提交修改 |
| `git push` | 推送到远程 |
| `git pull` | 拉取远程更新 |
| `git log --oneline` | 查看提交历史 |
| `git diff` | 查看未暂存的修改 |
| `gh auth status` | 查看 GitHub CLI 登录状态 |

---

## 一键推送脚本

创建快捷脚本 `push.sh`：

```bash
#!/bin/bash
cd /Users/mac/Documents/SCONE-DNA
git add .
git commit -m "${1:-Update}"
git push
echo "✅ 推送完成！"
```

使用方法：
```bash
chmod +x push.sh  # 首次需要添加执行权限
./push.sh "你的提交信息"
```

---

## 排除的文件夹

以下文件夹在 `.gitignore` 中配置，**不会被上传**：

- `1.0/` - 旧版本存档
- `experiment_results/` - 实验结果数据
- `__pycache__/` - Python 缓存

---

## 遇到问题？

### 查看认证状态
```bash
gh auth status
```

### 重新登录 GitHub CLI
```bash
gh auth login
# 选择 GitHub.com → HTTPS → Login with web browser
```

### 设置 Git 使用 GitHub CLI 认证
```bash
gh auth setup-git
```

### 认证失败 (403)
```bash
# 方法1：重新登录 gh
gh auth login

# 方法2：用 token 临时推送
git remote set-url origin https://cruan1991:TOKEN@github.com/cruan1991/SCONE-DNA.git
git push
git remote set-url origin https://github.com/cruan1991/SCONE-DNA.git
```

### 冲突
```bash
# 先拉取远程更新
git pull --rebase
# 解决冲突后
git push
```

### 撤销最后一次提交（未推送）
```bash
git reset --soft HEAD~1
```

---

## 项目结构

```
SCONE-DNA/
├── fsm_constraint.py          # FSM约束控制器
├── masked_arithmetic_codec.py # 带掩码算术编码
├── minimal_arithmetic_codec.py # 标准算术编码器
├── scone_fsm_arith.py         # 主API
├── README.md                  # 项目文档
├── requirements.txt           # 依赖说明
├── GITHUB_PUSH_GUIDE.md       # 本文档
├── scripts/                   # 实验脚本
│   ├── scone_experiments.py
│   ├── scone_ablation_experiment.py
│   ├── plot_metrics.py
│   ├── visualize_fsm_steering.py
│   └── ecc_simulation.py
├── 1.0/                       # [不上传] 旧版本
└── experiment_results/        # [不上传] 实验结果
```

---

*最后更新：2026-01-25*
*认证方式：GitHub CLI (gh auth)*
