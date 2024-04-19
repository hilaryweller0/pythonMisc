import numpy as np
import matplotlib.pyplot as plt


x = np.linspace(0,2,41)
dt = np.pi/4

fig = plt.figure()
timer = fig.canvas.new_timer(interval = 1000) # 1000 miliseconds delay
timer.add_callback(plt.close)

for i in range(3):
    y = np.sin(x*np.pi - i*dt)
    plt.plot(x, y)
    timer.start()
    plt.show()

