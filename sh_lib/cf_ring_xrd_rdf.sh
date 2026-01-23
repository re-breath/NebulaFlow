# 该脚本用来解析碳纤维材料，解析碳纤维的La,Lc,d002,5，6，7碳圆环的数量
# 该脚本主要使用了xrd.exe和cf_analyze程序，但是后面准备使用debyer程序来计算xrd
# 2025-12-25


# 该函数用来专门解析碳纤维材料，解析碳纤维的La,Lc,d002,5，6，7碳圆环的数量
analyze_cf(){

    #需要使用到的软件
    cf_exe="cf_analyze"
    xrd_exe="xrd"
    argfile="$1"

    which $cf_exe || (echo "Error: $cf_exe 未安装，请先安装 from https://github.com/kaushikljoshi/cf_analyze" && exit 1) 
    which $xrd_exe || (echo "Error: $xrd_exe 未安装，请先安装 from https://github.com/kaushikljoshi/cf_analyze" && exit 1)

    # 生成cf_analyze需要的参数文件
    cat << EOF > settings.txt
input_file   $argfile  #Name of the input xyz structure
table_flag      3    #1 = lammps connection table, 2 = reax connection table, 3 = build_distance_based table
ring_flag       1    #1 means identify carbon rings
void_flag       0    #1 means identify voids.
periodic_flag   1    #1 means system is periodic. Dims should be in input_file
grid_size       4.0  #Grid size for binning. Will be used only if table_flag is 1
bo_cut_off      0.3  # bo cut-off for identifying molecules for lammps/reax connection table
distance_cut_off 1.7 #This cut_off will be used to identify bonds based on interatomic distace. Default is 1.7
empty_bin_cut_off 1   #If any bin contains atoms less than or equal to this number, then it will be tagged as empty
EOF
    # 运行cf_analyze
    $cf_exe 

    # 生成xrd需要的参数文件
    cat << EOF > md.rc
UNIT 14 md.input old
UNIT 15 $argfile old
UNIT 16 RDF_L1_L2_na.data unknown
UNIT 17 SFACTOR_L1_L2_na.data unknown
EOF

    cat << EOF > md.input
  1811      - NB (Number of number of bins to calculate Q)
  3700      - Nbins (Number of number of bins to calculate g(r)
0.0045      - K_del (step for Q integration)
    66      - Rcell (Size of computational cell for g(r))	
     1      - NRCFLAG (!Nearest image convention flag (0 = Not used, 1 = used))
     0      - WFLAG (If window function is used (0 = Not used, 1 = sinc function is used, 2 = Hann function))	
1.5418      - LAMBDA (Wavelength to calculate in 2Theta, ANGSTROMS)
     1      - UNITFLAG (Units to calculate scattering angle (0 = Inverse angstroms, 1 = 2Theta))
     1      - AFACTORFLAG (Account for atomic scattering factor (0 = Calculate structure factor, 1 = Calculate intensity))
EOF
    cp $argfile opted_back

    # 提取第二行的 Lattice 信息并生成替换行
    new_line=$(awk 'NR==2 {
    lattice=$0
    sub(/.*Lattice="/,"",lattice)
    sub(/".*/,"",lattice)
    split(lattice,b," ")
    printf("0 %s 0 %s 0 %s 90 90 90", b[1], b[5], b[9])
    }' opted_back)

    # 用 sed 替换第二行
    sed "2c\\$new_line" opted_back > $argfile

    xrd
    cp ~/.rebreath/plot_library/plot_rdf_xrd.py ./
    python3 plot_rdf_xrd.py 
    cp ~/.rebreath/deal_data/xrdtreatment.py ./
    python3 xrdtreatment.py 

    # 打印5，6，7圆环原子数量 rings_size_distribution.dat
    awk '{print $6,$7,$8}' rings_size_distribution.dat
}
