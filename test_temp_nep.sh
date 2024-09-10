nepfile=nep.txt

nepfile_name=${nepfile%.*}


compute_elastic_moduli(){
    compute_lib=/home/dhk/.rebreath/compute_lib
    sed "s/nepfile/$1/g" $compute_lib/calorine_compute_elastic.py > calorine_compute_elastic.py
    python3 calorine_compute_elastic.py > elastic_calorine.txt
    rm -f calorine_compute_elastic.py 
}
#初始化
source_address=/home/dhk/rebreath/AlN/auto_series_work
dname=${nepfile_name}_auto_series_work
cp -r $source_address ./${dname}
cp $nepfile ./${dname}
echo "-------------------->已完成初始化"


echo "-------------------->已画出最终nep训练图"

cd ${dname}
workaddress=$PWD
cp $nepfile $workaddress/elastic/rs
cd $workaddress/elastic/rs
compute_elastic_moduli $nepfile 
echo "-------------------->rs的弹性模量"
cat $workaddress/elastic/rs/elastic_calorine.txt

cp $nepfile $workaddress/elastic/wz
cd $workaddress/elastic/wz
compute_elastic_moduli $nepfile
echo "-------------------->wz的弹性模量"
cat $workaddress/elastic/wz/elastic_calorine.txt
cd $workaddress
echo "-------------------->已完成弹性模量计算"
