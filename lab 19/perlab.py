import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 定義阻尼諧振子的微分方程系統
# 將二階微分方程 d^2x/dt^2 + beta dx/dt + x = 0 轉換為一階微分方程組：
# y[0] = x (位移)
# y[1] = v (速度)
# dy[0]/dt = y[1]
# dy[1]/dt = -2 * lambda_val * y[1] - y[0]
def damped_oscillator(t, y, lambda_val):
    dxdt = y[1]
    dvdt = -lambda_val * y[1] - y[0] + np.sin(t)
    return [dxdt, dvdt]

# 設定起始條件
x0 = 1.0  # 初始位移 x(0) = 1
v0 = 0.0  # 初始速度 v(0) = 0
initial_conditions = [x0, v0]

# 設定時間範圍
t_span = (0, 50)
t_eval = np.linspace(t_span[0], t_span[1], 500) # 在此時間範圍內產生500個點用於繪圖，使曲線平滑

# 讓使用者輸入三組不同的 lambda 數值
lambda_values = [0.5, 1.0, 1.5, 2.0, 2.5]

# 準備繪圖
plt.figure(figsize=(10, 6))

# 對每一組 lambda 數值求解並繪圖
for lambda_val in lambda_values:
    # 使用 solve_ivp 函數求解微分方程
    # method='RK45' 是一種常用的數值積分方法
    sol = solve_ivp(
        lambda t, y: damped_oscillator(t, y, lambda_val),
        t_span,
        initial_conditions,
        t_eval=t_eval,
        method='RK45'
    )
    # 繪製位移 x(t)
    plt.plot(sol.t, sol.y[0], label=f'beta = {lambda_val}')

# 設定圖表標題和座標軸標籤
plt.title('Block system with damping')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.tight_layout()
# plt.grid(True)
plt.xlim(0,40)
plt.legend()   # 顯示圖例 (各條曲線對應的 lambda 值)
plt.show()     # 顯示圖表