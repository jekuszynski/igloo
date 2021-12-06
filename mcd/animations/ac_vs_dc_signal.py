import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#setup x and y data
x = np.arange(0, 20, 0.1)
y1 = np.sin(x)
y2 = x * 0

#setup plot
fig, ax = plt.subplots(figsize=(8,6))

#define two lines for animation
line1, = ax.plot(x, y1, color = "red",lw=3,linestyle='-',zorder=5)
line2, = ax.plot(x, y2, color = "black",lw=3,linestyle='--',zorder=4)

#set plot limits
plt.xlim(0,20)
plt.ylim(-1.3,1.3)

# remove axes ticks
plt.xticks([],'')
plt.yticks([],'')

#function for updating using in animation
def update(num, x, y1, y2, line1, line2):
    line1.set_data(x[:num], y1[:num])
    line2.set_data(x[:num], y2[:num])
    return [line1,line2]

#add lines and text to plot
plt.plot([0,20],[1,1],color='black',lw=1,linestyle='--',zorder=3) #add y=1 AC reference line
plt.plot([0,20],[-1,-1],color='black',lw=1,linestyle='--',zorder=2) #add y=-1 AC reference line
plt.plot([0,20],[0,0],color='grey',lw=1,linestyle='--',zorder=1) #add y=0 DC reference line
plt.annotate(r'$V_{AC}$ Signal',(15.75,1.1)) #add V_AC Label
plt.annotate(r'$V_{DC}$ Signal',(15.75,0.1)) #add V_DC Label
ani = animation.FuncAnimation(fig, update, frames=len(x), fargs=[x, y1, y2, line1, line2], interval=20, blit=True, repeat=True, save_count=len(x))

#label x and y axes using r'' for raw text (LaTeX style)
ax.set_xlabel(r'Time ($t$) $\longrightarrow$') 
ax.set_ylabel(r'Signal ($V$)')

#set plot style
plt.style.use('seaborn-talk')

#show plot
plt.show()

#save animation
ani.save('ac_vs_dc_signal.mp4', dpi=300, savefig_kwargs={'transparent': True})
ani.save('ac_vs_dc_signal.gif', dpi=300, savefig_kwargs={'transparent': True})