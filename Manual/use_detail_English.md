## Detailed explanation of use (English version)

- [NEP related commands](#nep-related-commands)
- [Simplify complex operations](#simplify-complex-operations)
- [Format conversion](#format-conversion-related-commands)
- [GPUMD data processing](#gpumd-related-commands)
- [VASP calculation](#vasp-calculation)

## NEP related commands

```get_energy```
* Get all the energy in the xyz file
* The method used is get_energy $xyz_file, where $xyz_file is the xyz file path, and the output is the energy value of, if you do not enter the xyz file, it will default to get all the energy dump.xyz file configuration.


```get_Lattice```
* Get lattice parameters in xyz file
* The method used is get_Lattice $xyz_file, where $xyz_file is the xyz file path and the output is the lattice parameter.

```get_virial```
* Get the volume in the xyz file
* The method used is get_virial $xyz_file, where $xyz_file is the xyz file path, and the output is all the bit force information of the xyz file.


```get_configs_num```
* Get the number of configurations in the xyz file
* The method used is get_configs_num $xyz_file, where $xyz_file is the xyz file path, and the output is the number of configurations in the xyz file. The default xyz file is train.xyz.

```get_V```
* Get the volume in the xyz file
* The method used is get_V $xyz_file, where $xyz_file is the xyz file path and the output is the volume in the xyz file. The default xyz file is dump.xyz.

```screening_reasonable_forces```
* Filter the training dataset of nep and extract the reasonable configuration of the training dataset
* The method used is screening_reasonable_forces $xyz_file $min_force $max_force, where $xyz_file is the xyz file path, $min_force is the least reasonable force, $max_force is the maximum reasonable force, and the output is the filtered training dataset.

```screening_reasonable_energy```
* Screen the training dataset of the nep to extract the reasonable energy of the training dataset
* The method used is screening_reasonable_energy $xyz_file $min_energy $max_energy, where $xyz_file is the xyz file path, $min_energy is the least reasonable energy, $max_energy is the maximum reasonable energy, and the output is the filtered training dataset.

```screening_reasonable_virial```
* Filter the training dataset of the nep and extract the reasonable potential force of the training dataset
* The method used is screening_reasonable_virial $xyz_file $min_strain $max_strain, where $xyz_file is the xyz file path, $min_strain is the least reasonable bit force, $max_strain is the maximum reasonable bit force, and the output is the filtered training dataset.

```plot_nep```
* Draw the result graph of nep
* The method used is to use plot_nep directly in the directory of training nep, call python3 to draw, and output it as a picture file of.

## Simplify complex operations

```free_time_run```
* Perform tasks after monitoring idle GPUs
* The method used is the free_time_run 'command, where the command is the command that needs to be run, and when an idle gpu appears, the command will be run automatically.
  * For example: free_time_run'gpumd ', when an idle gpu appears, it will automatically run gpumd, which is equivalent to delaying gpumd's computation until a GPU is available.

```average_file```
* Calculate the average of multiple xyz files, using C++ implementation, and enter the average into the average.out file
* The method used is average_file $file1 $file2 $file3..., where $file1 $file2 $file3... is the file path and the output is the average in the file.
  * Use the example average_file thermo * to average all the thermo * files in the current folder and output as an average.out file.

```average_file_s```
Calculate the average across multiple xyz files, using a shell implementation
* The method used is average_file_s $file1 $file2 $file3..., where $file1 $file2 $file3... is the file path, and the output is the average of the corresponding data in the file, which is more flexible and will be automatically renamed.

```To replot```
* Paint function, which will draw the first two columns in the specified file
* The method used is replot $data_file, where $data_file is the data file path, the first column is x column data, the second column is y column data, call python3 to draw, and output as the picture file of.


```find_column_max```
* Find the maximum value of the specified column in the xyz file
* The method used is find_column_max $filename $column, where $filename is the xyz file path, $column is the column number, and the output is the maximum value for that column.


```find_column_abs_max```
* Find the absolute maximum value of the specified column in the xyz file
* The method used is find_column_abs_max $filename $column, where $filename is the xyz file path, $column is the column number, and the output is the absolute maximum value of that column.

## Format conversion related commands

```xyz_to_poscar```
* Convert xyz file to poscar file
* The method used is xyz_to_poscar $xyz_file, where $xyz_file is the xyz file path, and the output is the poscar file. (Note that the python library called for ovito needs to be prepared in advance)

```poscar_to_xyz```
* Convert poscar file to xyz file
* The method used is poscar_to_xyz $poscar_file, where $poscar_file is the poscar file path, and the output is the xyz file. (Note that the python library called for ovito needs to be prepared in advance)

```xyz_to_cssr```
* Convert xyz file to cssr file
* The method used is xyz_to_cssr $xyz_file $cssr_file, where $xyz_file is the xyz file path, $cssr_file is the cssr file path, and the output is the cssr file. (Note that the Python library called for ase needs to be prepared in advance)

```xyz_to_cif```
* Convert xyz file to cif file
* The method used is xyz_to_cif $xyz_file $cif_file, where $xyz_file is the xyz file path, $cif_file is the cif file path, and the output is the cif file. (Note that the Python library called for ase needs to be prepared in advance)

## gpumd related commands

```compute_elastic_moduli```
* Calculate elastic modulus using calorine
* The method used is compute_elastic_moduli $nep_file, where $nep_file is the nep file path, and the output is the elastic modulus of. (Note that the Python library called for calorine needs to be prepared in advance)


```plot_hnemd```
* Draw the result graph of hnemd
* The method used is to use plot_hnemd directly in the directory of training hnemd. The command will automatically determine the direction of hnemd's calculation of thermal conductivity, call python3 to draw, and output as a picture file of. (Note that the determination of this command is the name of the current directory. The name needs to have _x or _y or _z, otherwise an error will be reported)

```plot_mul_hnemd```
* Draw the result graph of multiple hnemds
* The method used is to use it directly in the directory where multiple hnemds are trained plot_mul_hnemd
  * For example, there are hnemd_0, hnemd_1, hnemd_2 three folders in the current directory, and the thermal conductivity map of the three folders will be automatically drawn after use. The command does not require other parameters. Note that it is automatically determined to determine whether there is _ [xyz] in the path of the current directory, so it is necessary to ensure the naming specification of the folder.

```deal_hnemd_data```
* Process the data completed by the hnemd calculation
* Similar to plot_mul_hnemd, but more perfect.
  * Venue requirements: The place where the function is used needs to have a lot of hnemd_ [0-9] + type folders


```cell_expansion```
* Use gpumd for cell expansion
* The method used is cell_expansion $nx $ny $nz, where $nx $ny $nz is the multiple of the expanded cell, and the output is the xyz file after the expanded cell. (Note that this command calls gpumd)
  * Since the cell expansion function of gpumd is mobilized, the command needs to have a nep.txt file in the default folder, and the model.xyz file is expanded by default
  * You can use xyzfile = test, nepfile = nep_test to modify the function's specified model with nep file


```plot_stress_strain_curve```
Automatically detect the strain of the axis, output data to the file, and draw (for drawing stress-strain curves).
* The method used is plot_stress_strain_curve
   * Requirements: You need to have _ {xyz} + file in the current directory name, which is the strain and stress file calculated by gpumd


```plot_mul_stress_strain_curve```
* Automatically detect multiple _ {xyz} + files in the current directory and automatically process their data to draw a picture
* The method used is plot_mul_stress_strain_curve
   * Suggestion: Name folders deform_x, deform_y, deform_z when uniaxial stretching


```relib```
* Administrative command, which is provided to users who have some knowledge of the NebulaFlow library, allowing them to quickly locate files in the library and place them in their current folder.
  * Option -v: View the file that will be found
     * For example: relib -v hp View all plot-related files in the NebulaFlow library
  * Option -m: Multi-level lookup
     * For example: relib -m hp pre multi-level lookup for files with both hp and pre in the filename of the NebulaFlow library


## VASP calculation

Before using the command module, you need to pay attention to learning the kspacing keyword to prevent inconsistencies in accuracy. All the following vasp calculation commands will discard kpoints and use kspacing uniformly.

Note: Remember to assign the value to kspacing before using the following command, for example: kspacing = 0.2

The following commands will automatically add kspacing to the INCAR file. If you do not assign a value in advance, 0.2 will be used by default.

```add_kspacing_to_incar```
* Add kspacing to incar file
* The method used is add_kspacing_to_incar kspacing, which defaults to 0.2.


```check_vasp_complete```
* Check if vasp completes the calculation
* The method used is check_vasp_complete, and if the calculation is not completed, vasp will be run.
  This command can replace the command submitted by VASP to some extent, preventing VASP from being run multiple times.
  * Note: This command defaults to using vasp submission, defaults to one graphics card or one core, if you need to use multiple cores, you need to add the number of cores after the command, for example: check_vasp_complete vasp 4

```run_all_vasp_job```
* Run vasp calculation (one core is used by default, if using cpu, you can add core_num = number of cores before using the command to specify the number of cpu cores used)
* The method used is run_all_vasp_job, which will look for the train- * folder in the current directory and run vasp.
  This command will automatically add kspacing to the incar file and check if vasp has completed its calculations.

```load_single_point_energy_dir```
Check all the POSCAR class files in the current directory (need to have POSCAR in the name) and organize them into a folder to facilitate the next calculation.
* The requirement is that the directory of the poscar file has incar and potcar below.
* Tip: You can use the run_all_vasp_job command for batch calculations immediately after using this command

```deal_outcar_to_train```
* Organize the OUTCAR file of the vasp calculation into a training dataset
* Find all current OUTCAR files and organize them into train.xyz
* Note that it is looking for all OUTCAR files, so you need to ensure that all OUTCAR files in the current directory are OUTCAR that can be calculated at a single point