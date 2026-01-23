#该脚本使用来安装NebulaFlow库

mkdir -p $HOME/.rebreath

rsync -av --update . $HOME/.rebreath/

if ! grep -q "WRITENEBULAFLOW2ENV" "$HOME/.bashrc"; then
    echo '
#rebreath的函数库
#------------------------------------------------

WRITENEBULAFLOW2ENV=1      
export PATH=$PATH:$HOME/.rebreath
source $HOME/.rebreath/rebreath-env-function

#------------------------------------------------
' >> "$HOME/.bashrc"
fi


echo "NebulaFlow库已经安装成功，下面请使用 source $HOME/.bashrc 命令使环境变量生效。"
echo "The NebulaFlow library has been successfully installed. Please use the source $HOME/.bashrc command to make the environment variables take effect."