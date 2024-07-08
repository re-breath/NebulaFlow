
![image](https://github.com/re-breath/NebulaFlow/blob/main/logo/NebulaFlow-logo.png)
## NebulaFlow  

## What is NebulaFlow?

* NebulaFlow是一款辅助分子动力学（MD）密度泛函理论（DFT）计算的命令行工具，内部设置了大量的命令来帮助完成计算与计算后的数据处理。

* NebulaFlow is a command-line tool that assists molecular dynamics (MD) first-principles (DFT) calculations. A large number of commands are set up internally to assist in the completion of calculations and post-computational data processing.

## What is the use of NebulaFlow?

* 帮助广大的科研工作者简化工作流程，以及降低处理数据的门槛，在我刚入门的时候，总是无法处理计算后的数据，面对github上很多库却不知道如何使用，使用NebulaFlow可以解决这个问题，只需要在命令行中敲出一个命令就能完成一系列的任务，无需再为不知道如何处理相关任务而烦恼。

* To help the majority of scientific researchers simplify their workflow and lower the threshold for processing data, when I first started, I was always unable to handle the calculated data. In the face of many libraries on GitHub, I did not know how to use them. Using NebulaFlow can solve this problem by simply typing a command on the command line to complete a series of tasks, without having to worry about not knowing how to handle related tasks.

* NebulaFlow的核心是使用shell来操控python，C++以及其他的语言，从而实现各种各样的功能，能极大的降低使用时的学习成本。

* The core of NebulaFlow is to use the shell to control Python, C++ and other languages, so as to achieve a variety of functions, which can greatly reduce the learning cost when using.

## NebulaFlow currently supports computing

* 目前支持GPUMD与VASP相关的计算，主要包括GPUMD构建训练集与训练NEP等的各种操作，以及VASP计算单点能等计算。
* At present, GPUMD and VASP-related calculations are supported, including various operations such as GPUMD building training datasets and training NEP, as well as VASP computing single point energy and other calculations.

* 本项目立志成为宇宙第一辅助数据处理的项目，以后会逐渐添加对更多计算的支持。
* This project aspires to become the first auxiliary data processing project in the universe, and will gradually add support for more calculations in the future.

## How to use NebulaFlow

NebulaFlow使用方法为在linux环境下输入各种命令即可，目的是使用一行简单的命令完成一系列复杂的任务，详细的使用方法如下：

[点击这里查看中文版使用手册](https://github.com/re-breath/NebulaFlow/blob/main/Manual/use_detail_Chinese.md)

NebulaFlow can be used to enter various commands in the Linux environment. The purpose is to complete a series of complex tasks with a single simple command. The detailed usage method is as follows:

[Click here to view the English version of the user manual](https://github.com/re-breath/NebulaFlow/blob/main/Manual/use_detail_English.md)

## Install NebulaFlow

* 使用NebulaFlowinstaller.sh脚本进行一键安装
* 安装方法： bash NebulaFlowinstaller.sh
  * 安装后使用 source $HOME/.bashrc，如果显示为 NebulaFlow library loaded O.<，并且无任何报错则代表安装成功。
  * 注意安装中可能会出现window转为linux后的格式问题而导致报错，可以使用dos2unix ~/.rebreath/*来尝试修复格式问题。
 
* Use NebulaFlowinstaller.sh script for one-click installation
* Installation method: bash NebulaFlowinstaller.sh
  * Use source $HOME/.bashrc after installation. If it shows NebulaFlow library loaded O. < and there is no error, the installation is successful.
  * Note that there may be a format problem after window is converted to linux during installation, resulting in an error. You can use dos2unix~/.rebreath/* to try to fix the format problem.
```shell
tar -zxvf NebulaFlow-1.0.tar.gz
cd NebulaFlow-1.0.tar.gz
bash NebulaFlowinstaller.sh
source ~/.bashrc
```
## 使用要求

* 要求linux中可以使用python3,内置命令默认使用python3。

* NebulaFlow使用shell调用了各种python库对文件进行处理，其很多命令无需任何python库也能实现，但是一些特殊的函数调用了部分库，只需要看起报错时候的提示需要哪些库即可。

* 可能用到的库有matplotlib,numpy,pylab,calorine,gpyumd,ase,gpumd-wizard等库，用到时候再进行下载即可

* Required to use python3 in linux, built-in commands use python3 by default.

NebulaFlow uses the shell to call various Python libraries to process files, and many of its commands can be implemented without any Python libraries. However, some special functions call some libraries, and you only need to see which libraries are required when the error is reported.

* Libraries that may be used include matplotlib, numpy, pylab, calorine, gpyumd, ase, gpump-wizard and other libraries. You can download them when you use them.

