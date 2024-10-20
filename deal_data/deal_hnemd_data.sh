 #!/bin/bash
 #本脚本的任务为：将进行多次HNEMD计算的结果进行平均处理，而后使用画出hnemd结果的图来简化流程

init_address=$PWD
mkdir  -p average_hnemd
mkdir  -p average_hnemd/kappa
mkdir  -p average_hnemd/shc

#整理各个文件到相应的地方方便管理，对齐颗粒度
function granularity(){
    local m=0
    file_column=0
    for i in $(find $init_address -maxdepth 1 -type d -regex  ".*hnemd[\._]?[0-9]+" | sort -t '_' -k 2 -n);do
        basei=$(basename $i)
        file_column=$(awk '{print NF}' $i/kappa.out)
        add_name=$(grep -oE '[0-9]+' <<< $basei | tail -n 1)
        if [ -z $add_name ]
        then
            add_name=$m
        fi
        m=$((m+1))
        cp $i/kappa.out $init_address/average_hnemd/kappa/kappa_${add_name}.out
        cp $i/thermo.out $init_address/average_hnemd/kappa/thermo_${add_name}.out
        cp $i/shc.out $init_address/average_hnemd/shc/shc_${add_name}.out
    done
}


#检查文件行和列是否一致
check_Rows_and_columns() {
    local files=("$@")
    local num_files=${#files[@]}
    local line_counts=()
    local column_count=0

    # Step 1 & 2: 检查文件存在性和行列一致性
    for file in "${files[@]}"; do
        if [ ! -f "$file" ]; then
            echo "文件 $file 不存在，请检查。" >&2
            exit 1
        fi
        
        local current_line_count=$(wc -l < "$file")
        line_counts+=($current_line_count)

        if [ $column_count -eq 0 ]; then
            column_count=$(awk '{print NF; exit}' "$file")
        elif [ $(awk '{print NF; exit}' "$file") -ne $column_count ]; then
            echo "文件 $file 的列数与之前的列数不匹配，请检查。" >&2
            exit 1
        fi
    done
    # 检查所有文件的行数是否一致
    for count in "${line_counts[@]}"; do
        if [ $count -ne ${line_counts[0]} ]; then
            echo "文件的行数不一致，请检查。" >&2
            exit 1
        fi
    done
}

# average_file(){                  
#     files_num="$#"
#     tmp_file='temp'
#     paste "$@" >$tmp_file
#     # 逐行读取处理
#     while IFS= read -r line
#     do
#         echo "$line" | awk -v files_num="$files_num" '{
#             sum=0;
#             for(i=1; i<=NF/files_num; i++){
#                 sum=0;
#                 for(j=0; j<files_num; j++){
#                     sum += $(i+j*NF/files_num)
#                 }
#                 printf "%.15f\t", sum/files_num;  
#             }
#             printf "\n";
#         }'
#     done < "$tmp_file" >average.out
#     rm -f "$tmp_file"
# }



granularity
cd $init_address/average_hnemd/kappa
check_Rows_and_columns

>kappa.out
average_file $(find ./ -type f -regex '.*kappa_[0-9]+\.out' |sort)
cp average.out kappa.out

# >thermo.out
# average_file $(find ./ -type f -regex '.*thermo_[0-9]+\.out' |sort)

