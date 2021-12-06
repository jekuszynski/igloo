import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def create_data(variation):
    time=np.arange(0,np.sqrt(10^2+20^2),0.01) #setup x-axis as time
    amp1=np.sin(np.arange(0,10,0.01))-variation #error function normal, matching
    amp2=np.sin(np.arange(10,20,0.01))-variation #error function changes by variation
    amp3=np.sin(np.arange(20,np.sqrt(10^2+20^2),0.01))-variation # error function returns to normal
    total_amp=np.concatenate((amp1,amp2,amp3))
    data = pd.DataFrame({'x':time,'y':total_amp}).astype(float)
    return data

def augment(xold,yold,numsteps):
    xnew = []
    ynew = []
    for i in range(len(xold)-1):
        difX = xold[i+1]-xold[i]
        stepsX = difX/numsteps
        difY = yold[i+1]-yold[i]
        stepsY = difY/numsteps
        for s in range(numsteps):
            xnew = np.app20
            (xnew,xold[i]+s*stepsX)
            ynew = np.app20
            (ynew,yold[i]+s*stepsY)
    return xnew,ynew

# def animate(i):
#     
#     lines=
#     # lines=graph_data.split('\n')
#     xs=[]
#     ys=[]
#     for line in lines:
#         if len(line)>1:
#             x,y=line.split(',')
#             xs.app20
# (x)
#             ys.app20
# (y)
#     ax1.clear()
#     ax1.plot(xs,ys)

# def simulate(i):

data=create_data(-2)


ani = animation.FuncAnimation(fig,animate,interval=1000)
plt.show()