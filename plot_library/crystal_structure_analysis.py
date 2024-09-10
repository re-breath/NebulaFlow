#本脚本的功能为对dump.xyz文件的每一帧进行晶型分析，分析其蜂窝状的原子的个数并画出图形进行可视化

import matplotlib.pyplot as plt
import matplotlib
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier
from ovito.data import *


# 设置matplotlib不显示图形界面
matplotlib.use('Agg')

# 函数：计算并绘制每个文件的石墨烯原子计数
def plot_graphene_counts(file_pattern, color, label):
    pipeline = import_file(file_pattern, multiple_frames=True)

    ptm_modifier = PolyhedralTemplateMatchingModifier()
    ptm_modifier.structures[PolyhedralTemplateMatchingModifier.Type.GRAPHENE].enabled = True
    pipeline.modifiers.append(ptm_modifier)

    graphene_counts = []
    frames = []

    for frame_index in range(pipeline.source.num_frames):
        data = pipeline.compute(frame_index)
        frames.append(frame_index)
        graphene_counts.append(data.attributes['PolyhedralTemplateMatching.counts.GRAPHENE'])

    plt.plot(frames, graphene_counts, 'o-', color=color, label=label,markersize=2)

# 函数：设置图表并保存为PNG文件
def setup_and_save_plot(file_list):
    plt.figure(figsize=(10, 5))
    for file_data in file_list:
        plot_graphene_counts(file_data['pattern'], file_data['color'], file_data['label'])

    plt.title('Graphene Atom Count Over Time')
    plt.xlabel('Time/ns$^{-1}$')
    plt.ylabel('Count of Graphene Atoms')
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.savefig("graphene_atom_count.png")

# 文件列表和颜色设置
file_list = [
    {'pattern': 'dump.xyz', 'color': 'blue', 'label': 'down-500k-1ns'},
    #{'pattern': 'dump_250k_1ns.xyz', 'color': 'green', 'label': 'down-250k-1ns'},
    # 更多文件可以添加到列表
]

# 调用函数
setup_and_save_plot(file_list)

