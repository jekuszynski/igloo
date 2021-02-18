import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
 
# References
# https://towardsdatascience.com/animations-with-matplotlib-d96375c5442c
# https://riptutorial.com/matplotlib/example/23558/basic-animation-with-funcanimation
 
def func(t, line):
    t = np.arange(0,t,0.1)
    y = np.sin(t)
    line.set_data(t, y)
    return line

fig = plt.figure()
ax = plt.axes(xlim=(0, 20), ylim=(-1.2, 1.22))
redDots = plt.plot([], [], 'ro')
line = plt.plot([], [], lw=2)

# Creating the Animation object
line_ani = animation.FuncAnimation(fig, func, frames=np.arange(1,20,0.1), fargs=(line), interval=20, blit=False)
#line_ani.save(r'Animation.mp4')

plt.show()

line_ani.save('live_plot_sine_wave.gif')