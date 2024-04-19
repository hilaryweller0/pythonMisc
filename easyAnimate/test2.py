import matplotlib.pyplot as plt


fig = plt.figure()
timer = fig.canvas.new_timer(interval = 1000) # 1000 miliseconds delay
timer.add_callback(plt.close)

plt.plot([1,2,3,4])
plt.ylabel('some numbers')

timer.start()
plt.show()
print("Am doing something else")
