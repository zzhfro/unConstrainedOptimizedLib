
import numpy as np
import matplotlib.pyplot as plt

# 定义二维函数
def f(x, y):

  result=0.01*x*x+y*y 
  return result

# 生成x和y的网格
x = np.linspace(-10, 10, 1000)
y = np.linspace(-10, 10, 1000)

x, y = np.meshgrid(x, y)

# 计算二维函数的值
z = f(x, y)

# 绘制二维函数图形
levels = np.linspace(0, 50, 50)  # 设置50个等高线级别
plt.contourf(x, y, z, levels=levels)
plt.colorbar(label='Function Value')
plt.title('2D Function Plot')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show()
