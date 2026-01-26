# GitHub 推送指南

## 仓库信息

| 项目 | 值 |
|------|-----|
| **仓库地址** | https://github.com/cruan1991/SCONE-DNA |
| **用户名** | cruan1991 |
| **本地路径** | /Users/mac/Documents/SCONE-DNA |

---

## 日常更新流程

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

## 首次推送认证

首次推送时会要求输入凭证：

```
Username: cruan1991
Password: [你的 Personal Access Token]
```

⚠️ **注意**：Password 填的是 **Personal Access Token**，不是 GitHub 密码！

### 获取 Personal Access Token

1. 登录 GitHub
2. 点击右上角头像 → **Settings**
3. 左侧菜单最下方 → **Developer settings**
4. **Personal access tokens** → **Tokens (classic)**
5. **Generate new token**
6. 勾选 `repo` 权限
7. 生成并复制 token

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

---

## 一键推送脚本

可以创建一个快捷脚本 `push.sh`：

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

### 认证失败 (403)
```bash
# 重新设置 remote URL（带 token）
git remote set-url origin https://cruan1991:TOKEN@github.com/cruan1991/SCONE-DNA.git
git push
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
