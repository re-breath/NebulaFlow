#该脚本使用来安装NebulaFlow库

mkdir -p $HOME/.rebreath

rsync -av --update . $HOME/.rebreath/

echo 'export PATH=$PATH:$HOME/.rebreath' >> $HOME/.bashrc

echo '
#rebreath的函数库
#------------------------------------------------

source $HOME/.rebreath/rebreath-env-function

#------------------------------------------------
' >> $HOME/.bashrc

echo "NebulaFlow库已经安装成功，下面请使用 source $HOME/.bashrc 命令使环境变量生效。"
echo "The NebulaFlow library has been successfully installed. Please use the source $HOME/.bashrc command to make the environment variables take effect."
