## 使用详解（中文版）
- [NEP相关命令](#nep相关命令)
- [简化复杂操作](#简化复杂操作)
- [格式转换相关命令](#格式转化相关命令)
- [gpumd数据处理](#gpumd相关命令)
- [VASP计算相关命令](#VASP计算)

## NEP相关命令

```get_energy```
* 获得xyz文件中的所有的能量
* 使用的方法为get_energy $xyz_file，其中$xyz_file为xyz文件路径，输出为的能量值,如果如果不输入xyz文件则会默认为获取dump.xyz文件构型的所有的能量。


```get_Lattice```
* 获得xyz文件中的晶格参数
* 使用的方法为get_Lattice $xyz_file，其中$xyz_file为xyz文件路径，输出为的晶格参数。

```get_virial```
* 获得xyz文件中的体积
* 使用的方法为get_virial $xyz_file，其中$xyz_file为xyz文件路径，输出为的xyz文件所有的位力信息。


```get_configs_num```
* 获得xyz文件中的构型数量
* 使用的方法为get_configs_num $xyz_file，其中$xyz_file为xyz文件路径，输出为的xyz文件中的构型数量。默认xyz文件为train.xyz。

```get_V```
* 获得xyz文件中的体积
* 使用的方法为get_V $xyz_file，其中$xyz_file为xyz文件路径，输出为的xyz文件中的体积。默认xyz文件为dump.xyz。

```screening_reasonable_forces```
* 筛选nep的训练集，将训练集的合理的构型提取出
* 使用的方法为 screening_reasonable_forces $xyz_file $min_force $max_force，其中$xyz_file为xyz文件路径，$min_force为最小合理的力，$max_force为最大合理的力，输出为筛选后的训练集。

```screening_reasonable_energy```
* 筛选nep的训练集，将训练集的合理的能量提取出
* 使用的方法为 screening_reasonable_energy $xyz_file $min_energy $max_energy，其中$xyz_file为xyz文件路径，$min_energy为最小合理的能量，$max_energy为最大合理的能量，输出为筛选后的训练集。

```screening_reasonable_virial```
* 筛选nep的训练集，将训练集的合理的位力提取出
* 使用的方法为 screening_reasonable_virial $xyz_file $min_strain $max_strain，其中$xyz_file为xyz文件路径，$min_strain为最小合理的位力，$max_strain为最大合理的位力，输出为筛选后的训练集。

```plot_nep```
* 画出nep的结果图
* 使用的方法为在训练nep的目录下直接使用 plot_nep，调用python3画图，输出为的图片文件。

## 简化复杂操作

```free_time_run```
* 监测到空闲的gpu后进行任务
* 使用的方法为 free_time_run '命令，其中命令为需要运行的命令，当空闲的gpu出现时，会自动运行命令。
  * 例如：free_time_run 'gpumd'，当空闲的gpu出现时，会自动运行gpumd,相当于延后gpumd的计算直到有可用GPU。

```average_file```
* 计算多个xyz文件中的平均值,使用C++实现，平均后输入到average.out文件中
* 使用的方法为 average_file $file1 $file2 $file3...，其中$file1 $file2 $file3...为文件路径，输出为的文件中的平均值。
  * 使用实例average_file thermo* 将当前文件夹所有的thermo*文件进行平均，输出为average.out文件。

```average_file_s```
* 计算多个xyz文件中的平均值，使用shell实现
* 使用的方法为 average_file_s $file1 $file2 $file3...，其中$file1 $file2 $file3...为文件路径，输出为的文件中对应数据的平均值，更加灵活一些，会自动改名。

```replot```
* 画图函数，将指定文件中的前两列画图
* 使用的方法为 replot $data_file，其中$data_file为数据文件路径，第一列为x列数据，第二列为y列数据，调用python3画图，输出为的图片文件。


```find_column_max```
* 找到xyz文件中指定列的最大值
* 使用的方法为 find_column_max $filename $column，其中$filename为xyz文件路径，$column为列号，输出为该列的最大值。


```find_column_abs_max```
* 找到xyz文件中指定列的绝对值最大值
* 使用的方法为 find_column_abs_max $filename $column，其中$filename为xyz文件路径，$column为列号，输出为该列的绝对值最大值。

## 格式转化相关命令

```xyz_to_poscar```
* 将xyz文件转换为poscar文件
* 使用的方法为 xyz_to_poscar $xyz_file，其中$xyz_file为xyz文件路径，输出为的poscar文件。（注意调用的为ovito的python库，需要提前备好）

```poscar_to_xyz```
* 将poscar文件转换为xyz文件
* 使用的方法为 poscar_to_xyz $poscar_file，其中$poscar_file为poscar文件路径，输出为的xyz文件。（注意调用的为ovito的python库，需要提前备好）

```xyz_to_cssr```
* 将xyz文件转换为cssr文件
* 使用的方法为 xyz_to_cssr $xyz_file $cssr_file，其中$xyz_file为xyz文件路径，$cssr_file为cssr文件路径，输出为的cssr文件。（注意调用的为ase的python库，需要提前备好）

```xyz_to_cif```
* 将xyz文件转换为cif文件
* 使用的方法为 xyz_to_cif $xyz_file $cif_file，其中$xyz_file为xyz文件路径，$cif_file为cif文件路径，输出为的cif文件。（注意调用的为ase的python库，需要提前备好）

## gpumd相关命令

```compute_elastic_moduli```
* 使用calorine计算弹性模量
* 使用的方法为 compute_elastic_moduli $nep_file，其中$nep_file为nep文件路径，输出为的弹性模量。（注意调用的为calorine的python库，需要提前备好）


```plot_hnemd```
* 画出hnemd的结果图
* 使用的方法为在训练hnemd的目录下直接使用 plot_hnemd，该命令将会自动判断hnemd计算热导率的方向，调用python3画图，输出为的图片文件。（注意该命令的判定为当前目录的名字，名字中需要有_x或_y或_z，否则会报错）

```plot_mul_hnemd```
* 画出多个hnemd的结果图
* 使用的方法为在训练多个hnemd的目录下直接使用 plot_mul_hnemd
  * 例如当前的目录下存在hnemd_0,hnemd_1,hnemd_2三个文件夹，使用后则会自动画出三个文件夹的热导率图,该命令无需其他参数，注意自动判定为判定当前目录的路径中有没有_[xyz],因此需要保证文件夹的命名规范。

```deal_hnemd_data```
* 处理hnemd计算完成的数据
* 使用起来和plot_mul_hnemd类似，但更加完善。
  * 场地要求：使用该函数的地方需要有很多的hnemd_[0-9]+类型的文件夹


```cell_expansion```
* 使用gpumd进行扩胞
* 使用的方法为 cell_expansion $nx $ny $nz，其中$nx $ny $nz为扩胞的倍数，输出为的扩胞后的xyz文件。（注意该命令调用了gpumd）
  * 由于调动了gpumd的扩胞功能，因此该命令默认文件夹下需要有nep.txt文件，默认对model.xyz文件进行扩胞
  * 可以使用使用xyzfile=test , nepfile=nep_test.txt 修改该函数的指定的模型与nep文件


```plot_stress_strain_curve```
* 自动检测应变的轴，输出数据到文件中，并进行画图（用于画出应力应变曲线）
* 使用的方法为 plot_stress_strain_curve
   * 要求：需要在当前目录名字有_{xyz}+文件，该文件为gpumd计算的应变与应力文件


```plot_mul_stress_strain_curve```
* 自动检测当前目录下的多个_{xyz}+文件并自动将其数据进行处理画出图来
* 使用的方法为 plot_mul_stress_strain_curve
   * 建议：单轴拉伸时文件夹命名为deform_x,deform_y,deform_z


```relib```
* 管理命令,该命令提供给对NebulaFlow库有一定的了解的用户使用，可以快速查找库中的文件，并讲该文件放到当前文件夹下。
  * 选项  -v ：查看将会查找的文件
     * 例如：relib -v  hp 查看NebulaFlow库所有plot相关的文件
  * 选项  -m ：多级查找
     * 例如：relib -m hp pre  多级查找NebulaFlow库中文件名中既有hp又有pre的文件


## VASP计算

使用命令模块前需要注意学习kspacing关键字，防止出现精度不一致的情况，下面所有的vasp计算命令将会摒弃kpoints，统一使用kspacing。

注意：使用以下命令前记得先赋值给kspacing，例如：kspacing=0.2

下面的命令有些会自动的给INCAR文件添加kspacing，如果不提前赋值，将会默认使用0.2。

```add_kspacing_to_incar```
* 添加kspacing到incar文件中
* 使用的方法为 add_kspacing_to_incar kspacing，其中默认为0.2。


```check_vasp_complete```
* 检查vasp是否完成计算
* 使用的方法为 check_vasp_complete，如果没有完成计算，则会运行vasp。
  * 该命令可以一定程度上替代VASP提交的命令，防止多次运行VASP。
  * 注意：该命令默认为使用vasp提交，默认一个显卡或一个核，如果需要使用多个核，则需要在命令后面加上核数，例如：check_vasp_complete vasp 4

```run_all_vasp_job```
* 运行vasp计算（默认使用一个核数，如果使用cpu，可以在使用命令前添加 core_num=核数  来指定使用的cpu核数）
* 使用的方法为 run_all_vasp_job，该命令会在当前目录下查找train-*文件夹，并运行vasp。
  * 该命令会自动添加kspacing到incar文件中，并检查vasp是否完成计算。

```load_single_point_energy_dir```
* 检查当前目录下所有的POSCAR类（名字中需要有POSCAR）的文件，将其整理到一个文件夹中，方便进行下一步的计算。
* 要求是poscar文件的目录下面有incar与potcar
* 提示：使用该命令后可以立即使用run_all_vasp_job命令进行批量计算

```deal_outcar_to_train```
* 将vasp计算的OUTCAR文件整理成训练集
* 寻找当前所有的OUTCAR文件，将其整理成train.xyz
* 注意其寻找的是所有的OUTCAR文件，因此需要确保当前目录下所有OUTCAR文件都为单点能计算的OUTCAR


