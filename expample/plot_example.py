import matplotlib.pyplot as plt
import numpy as np

# x = np.loadtxt('data.txt')[:,0]
# y = np.loadtxt('data.txt')[:,1]

x = [0, 1, 2, 3, 4, 5]
y = [0, 1, 4, 9, 16, 25]

# 创建条形图
plt.figure(figsize=(10, 6))  # 设置图形大小
bars = plt.bar(x, y, color='skyblue', width=0.5)  # 设置颜色和条形宽度

# 添加标题和轴标签
plt.title('Bar Chart Example', fontsize=16, fontweight='bold')
plt.xlabel('X-axis Label', fontsize=12)
plt.ylabel('Y-axis Label', fontsize=12)

# 添加网格线
plt.grid(axis='y', linestyle='--', alpha=0.7)

# 添加数据标签
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval, round(yval, 2), ha='center', va='bottom', fontsize=10)

# 保存图表
plt.tight_layout()  # 自动调整布局
plt.savefig('test.png', dpi=300)  # 设置保存的图片分辨率

# 显示图表（可选）
# plt.show()