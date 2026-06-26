#!/usr/bin/env bash
# ======================================================================
# File:         NebulaFlowinstaller.sh
# Project:      NebulaFlow
# Description:  NebulaFlow一键安装脚本 / One-click installer for NebulaFlow.
#               Copies all files to ~/.rebreath/ and adds source to ~/.bashrc.
# Author:       rebreath
# ======================================================================

set -e

# -------------------------------------------------------------------------
# Step 1: Create library directory and sync files
# Step 1: 创建库目录并同步文件
# -------------------------------------------------------------------------
mkdir -p $HOME/.rebreath
rsync -av --update . $HOME/.rebreath/

# -------------------------------------------------------------------------
# Step 2: Add to ~/.bashrc if not already present
# Step 2: 将环境变量写入 ~/.bashrc（首次安装时）
# -------------------------------------------------------------------------
if ! grep -q "WRITENEBULAFLOW2ENV" "$HOME/.bashrc"; then
    echo '
# NebulaFlow function library by rebreath
# ----------------------------------------------------------------
WRITENEBULAFLOW2ENV=1
export PATH=$PATH:$HOME/.rebreath
export BASH_ENV="$HOME/.rebreath/rebreath-env-function"
source $HOME/.rebreath/rebreath-env-function
# ----------------------------------------------------------------
' >> "$HOME/.bashrc"
fi

# -------------------------------------------------------------------------
# Step 3: Make gpuq executable
# Step 3: 设置 gpuq 为可执行
# -------------------------------------------------------------------------
chmod +x $HOME/.rebreath/gpuq

# -------------------------------------------------------------------------
# Step 4: Done — 提示用户加载环境
# Step 4: Remind user to source ~/.bashrc
# -------------------------------------------------------------------------
echo ""
echo "NebulaFlow library installed successfully!"
echo "Please run: source ~/.bashrc"
echo "If you see 'NebulaFlow library loaded O.<', the installation is complete."
echo ""
echo "NebulaFlow库已经安装成功，下面请使用 source $HOME/.bashrc 命令使环境变量生效。"
echo "看到 'NebulaFlow library loaded O.<' 代表安装成功。"
