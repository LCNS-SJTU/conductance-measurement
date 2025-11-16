# !nrnivmodl /home/wzl/LFPy/project/conductance_measurement/realistic_neuron/hippocampal_pyramidal_neuron/mod
from neuron import h
import plotly.graph_objects as go
import matplotlib.colors as mcolors
from neuron.units import ms, mV
import plotly.io as pio
from matplotlib import pyplot, cm
import numpy as np
import matplotlib.pyplot as plt
# %matplotlib inline
import plotly
from neuron import clear_gui_callback
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from numpy.linalg import pinv, eig
import numpy as np
from multiprocessing import Pool
font2={'family':'Times New Roman',
'weight':'bold',
'size': 10}
# 自定义刻度标签显示格式

from matplotlib.ticker import FuncFormatter
from neuron import h, gui
import os

# 要更改的目标文档的路径
document_path = "/home/wzl/LFPy/project/conductance_measurement/realistic_neuron/hippocampal_pyramidal_neuron/conductance"

# 更改当前工作目录
os.chdir(document_path)

# 加载基础的 NEURON GUI 库
h.load_file("nrngui.hoc")


# 依次打开各个 HOC 文件
h.load_file("n128.hoc")               # 几何文件
h.load_file("axon_sections.hoc")      # 轴突部分
h.load_file("basal_dendrite.hoc")     # 基础树突
h.load_file("apical_dendrite.hoc")    # 顶端树突
h.load_file("apical_trunk.hoc")       # 顶端树干
h.load_file("radiatum.hoc")           # 放射区
h.load_file("init.hoc")               # 初始化设置
h.load_file("addgraph.hoc")           # 添加图形显示  


from neuron import h
from matplotlib import pyplot as plt

h.soma.diam = 1.
h.scale_bar.diam = 0.1
# h.dend[198].diam = 5.
# 绘图设置
ps = h.PlotShape(False)  # 创建PlotShape对象
ps.show(0)  # 不自动显示NEURON的窗口

# 使用Matplotlib进行绘图
ax = ps.plot(plt)

ax.view_init(elev = 90, azim = 105.6)
ax.grid(False)
# 关闭坐标轴
ax.axis('off')
plt.savefig("model_visualization.jpg", format='jpg', dpi=300)
plt.show()
