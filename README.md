
![image](https://github.com/re-breath/NebulaFlow/blob/main/logo/NebulaFlow-logo.png)
## NebulaFlow  

## What is NebulaFlow?

* NebulaFlow是一款辅助分子动力学（MD）密度泛函理论（DFT）计算的命令行工具，内部设置了大量的命令来帮助完成计算与计算后的数据处理。

* NebulaFlow is a command-line tool that assists molecular dynamics (MD) first-principles (DFT) calculations. A large number of commands are set up internally to assist in the completion of calculations and post-computational data processing.

## NebulaFlow有什么用

* 帮助广大的科研工作者简化工作流程，以及降低处理数据的门槛，在我刚入门的时候，总是无法处理计算后的数据，面对github上很多库却不知道如何使用，使用NebulaFlow可以解决这个问题，只需要在命令行中敲出一个命令就能完成一系列的任务，无需再为不知道如何处理相关任务而烦恼。

* To help the majority of scientific researchers simplify their workflow and lower the threshold for processing data, when I first started, I was always unable to handle the calculated data. In the face of many libraries on GitHub, I did not know how to use them. Using NebulaFlow can solve this problem by simply typing a command on the command line to complete a series of tasks, without having to worry about not knowing how to handle related tasks.

* NebulaFlow的核心是使用shell来操控python，C++以及其他的语言，从而实现各种各样的功能，能极大的降低使用时的学习成本。

* The core of NebulaFlow is to use the shell to control Python, C++ and other languages, so as to achieve a variety of functions, which can greatly reduce the learning cost when using.

## NebulaFlow目前支持的计算

* 目前支持GPUMD与VASP相关的计算，主要包括GPUMD构建训练集与训练NEP等的各种操作，以及VASP计算单点能等计算。
* At present, GPUMD and VASP-related calculations are supported, including various operations such as GPUMD building training datasets and training NEP, as well as VASP computing single point energy and other calculations.

* 本项目立志成为宇宙第一辅助数据处理的项目，以后会逐渐添加对更多计算的支持。
* This project aspires to become the first auxiliary data processing project in the universe, and will gradually add support for more calculations in the future.

## NebulaFlow 使用方法

* 筛选NEP训练集
  * ```screening_reasonable_forces```
    * 筛选nep的训练集，将训练集的合理的构型提取出，使用的方法为 screening_reasonable_forces xyzfile min_force max_force
    * Filter nep's training dataset, extract the reasonable configuration of the training dataset, using the method of screening_reasonable_forces xyzfile min_force max_force
      
  * ```screening_reasonable_energy```
    * 筛选nep的训练集，将训练集的合理的能量提取出,使用的方法为 screening_reasonable_energy xyzfile min_energy max_energy


  * ```screening_reasonable_virial```
    * 筛选nep的训练集，将训练集的合理的位力提取出,使用的方法为 screening_reasonable_virial xyzfile min_strain max_strain
    * Filter the training dataset of nep, and extract the reasonable potential force of the training dataset, using the method of screening_reasonable_virial xyzfile min_strain max_strain
      
