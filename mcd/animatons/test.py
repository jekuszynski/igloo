import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy import signal

plt.style.use('seaborn-colorblind') #set plot style

fig = plt.figure() #build figure
ax = plt.axes(xlim=(0, 12), ylim=(-1.2, .2)) #set plot axes
line, = ax.plot([], [], lw=2) #prepare the animated line

def init():
    line.set_data([], [])
    return line,

def animate(i): #function for animation process
    x = np.linspace(-2, 14, 800)
    # y = np.sin(2 * np.pi * (x - 0.01 * i))
    # y1 = [(np.sin(4 * np.pi * (x - 0.01 * i)))/8 for x in x if x <= 4] #default sin(x) centered at y=0
    y1 = [(0*x*i) for x in x if x <= 4] #line at y=0
    y2 = [((np.exp(-(x - 4.0)) - 1) + np.sin(4 * np.pi * (x - 0.01 * i))/8) for x in x if x > 4 and x <= 8] #exp function to show decrease in signal
    y3 = [((-np.exp(-(x - 8.0))) + np.sin(4 * np.pi * (x - 0.01 * i))/8) for x in x if x > 8] #exp function to show recovery of signal after correction
    y_convolution = np.concatenate((y1,y2,y3))
    y_smooth = signal.savgol_filter(y_convolution,23,3)
    # line.set_data(x, y)
    line.set_data(x, y_smooth)
    return line,

plt.plot([0,12],[0,0],color='k',lw=3,linestyle='dashed') #add y=0 reference line
plt.plot([0,12],[-0.02,-0.02],color='k',lw=1,linestyle='dashed') #add y=-0.1 reference line (error signal)
plt.plot([4,4],[-2,2],color='k',lw=2,linestyle='dashed') #add x=4 reference line
plt.plot([8,8],[-2,2],color='k',lw=2,linestyle='dashed') #add x=8 reference line

#try to adjust frames and interval to make as smooth as possible
anim = FuncAnimation(fig, animate, init_func=init, frames=100, interval=10, blit=True) #animate plot

plt.show()
anim.save('test_sine_wave.gif',fps=60,dpi=300,bitrate=1e6)