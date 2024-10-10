
generate_band_plot() {
    #函数说明：该函数用来指定DFPT的band.yaml文件拿来画图
    # 检查是否传入了 band.yaml 文件路径
    if [ $# -ne 1 ]; then
        echo "用法: generate_band_plot <band.yaml文件路径>"
        return 1
    fi

    band_file=$1

    # 检查文件是否存在
    if [ ! -f "$band_file" ]; then
        echo "文件不存在: $band_file"
        return 1
    fi

    # 调用Python脚本处理band.yaml并使用matplotlib绘图
    python3 <<EOF
import yaml
import numpy as np
import matplotlib.pyplot as plt

# 读取YAML文件
with open("$band_file", 'r') as file:
    data = yaml.safe_load(file)

# 提取q点的路径
q_points = [entry['distance'] for entry in data['phonon']]

# 提取所有频率
frequencies = []
for entry in data['phonon']:
    freqs = [band['frequency'] for band in entry['band']]
    frequencies.append(freqs)

# 转置列表以将解包频率与其对应的q点对齐
frequencies = list(map(list, zip(*frequencies)))

# 绘制每个频率
for freq in frequencies:
    plt.plot(q_points, freq)

# 设置坐标轴标签和标题
plt.xlabel('Wave Vector')
plt.ylabel('Frequency (THz)')
plt.title('Phonon Band Structure')
plt.savefig('phonon_band_structure.png')
plt.show()
EOF

    echo "声子谱图生成完毕。"
}

