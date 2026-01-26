# GitHub Push Guide

## Repository Information

| Item | Value |
|------|-------|
| **Repository URL** | https://github.com/cruan1991/SCONE-DNA |
| **Username** | cruan1991 |
| **Local Path** | /Users/mac/Documents/SCONE-DNA |
| **Authentication** | GitHub CLI (gh) ✅ Configured |

---

## Daily Update Workflow (One-Click Push)

```bash
cd /Users/mac/Documents/SCONE-DNA
git add .
git commit -m "Your update description"
git push
```

**No password or token required!** GitHub CLI is configured for automatic authentication.

---

## Detailed Steps

### 1. Navigate to Project Directory
```bash
cd /Users/mac/Documents/SCONE-DNA
```

### 2. Check Modification Status
```bash
git status
```

### 3. Add Modified Files
```bash
# Add all modifications
git add .

# Or add specific files
git add filename.py
```

### 4. Commit Changes
```bash
git commit -m "Brief description of your changes"
```

### 5. Push to GitHub
```bash
git push
```

---

## Command Quick Reference

| Command | Description |
|---------|-------------|
| `git status` | Check current status |
| `git add .` | Add all modifications |
| `git commit -m "msg"` | Commit changes |
| `git push` | Push to remote |
| `git pull` | Pull remote updates |
| `git log --oneline` | View commit history |
| `git diff` | View unstaged changes |
| `gh auth status` | Check GitHub CLI login status |

---

## One-Click Push Script

Create a shortcut script `push.sh`:

```bash
#!/bin/bash
cd /Users/mac/Documents/SCONE-DNA
git add .
git commit -m "${1:-Update}"
git push
echo "✅ Push complete!"
```

Usage:
```bash
chmod +x push.sh  # First time: add execute permission
./push.sh "Your commit message"
```

---

## Excluded Folders

The following folders are configured in `.gitignore` and **will not be uploaded**:

- `1.0/` - Old version archive
- `experiment_results/` - Experiment result data
- `__pycache__/` - Python cache

---

## Troubleshooting

### Check Authentication Status
```bash
gh auth status
```

### Re-login to GitHub CLI
```bash
gh auth login
# Select GitHub.com → HTTPS → Login with web browser
```

### Set Git to Use GitHub CLI Authentication
```bash
gh auth setup-git
```

### Authentication Failed (403)
```bash
# Method 1: Re-login to gh
gh auth login

# Method 2: Temporary push with token
git remote set-url origin https://cruan1991:TOKEN@github.com/cruan1991/SCONE-DNA.git
git push
git remote set-url origin https://github.com/cruan1991/SCONE-DNA.git
```

### Conflicts
```bash
# Pull remote updates first
git pull --rebase
# After resolving conflicts
git push
```

### Undo Last Commit (Not Pushed)
```bash
git reset --soft HEAD~1
```

---

## Project Structure

```
SCONE-DNA/
├── fsm_constraint.py          # FSM constraint controller
├── masked_arithmetic_codec.py # Masked arithmetic coding
├── minimal_arithmetic_codec.py # Standard arithmetic encoder
├── scone_fsm_arith.py         # Main API
├── README.md                  # Project documentation
├── requirements.txt           # Dependency description
├── GITHUB_PUSH_GUIDE.md       # This document
├── scripts/                   # Experiment scripts
│   ├── scone_experiments.py
│   ├── scone_ablation_experiment.py
│   ├── plot_metrics.py
│   ├── visualize_fsm_steering.py
│   └── ecc_simulation.py
├── 1.0/                       # [Not uploaded] Old version
└── experiment_results/        # [Not uploaded] Experiment results
```

---

*Last updated: 2026-01-26*
*Authentication: GitHub CLI (gh auth)*
