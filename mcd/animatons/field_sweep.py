import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import ticker

# Use matplotlib ggplot stylesheet if available
try:
    plt.style.use('ggplot')
except:
    pass

# Set which type of animation will be plotted. One of:
# line, pcolor, scatter, contour, quiver, labels
animation_type = 'line'

# ----------------------------------------------------------------------------
# Create data to plot. F is 2D array.

# Create a two-dimensional array of data: F(x, t)
xlim_end=np.pi*2
x = np.linspace(-0.5, xlim_end, 91)
b = np.linspace(0, 10, 60)
X2, B2 = np.meshgrid(x, b)
sinB2 = np.sin(B2*np.pi)
F = B2*np.sin(X2+0.14159)

# ----------------------------------------------------------------------------
# Set up the figure and axis
fig, ax = plt.subplots(figsize=(4, 2))

# ----------------------------------------------------------------------------
ax.set(xlim=(0-0.14159, xlim_end-0.14159), ylim=(-11, 11),xlabel=(r"$\omega_p$"), ylabel=(r'$\Delta A$ (x $10^{-3}$)'))
ax.xaxis.set_major_locator(ticker.NullLocator())

line = ax.plot(x, F[0, :], color='black', lw=2)[0]
plt.plot([3,3],[-15,15],color='black',linestyle='--',linewidth=2)

def animate(i):
    line.set_ydata(F[i, :])

plt.tight_layout()

# Save the animation
anim = FuncAnimation(fig, animate, interval=30, frames=len(b)-1, repeat=True)
fig.show()
anim.save('field_sweep.gif', dpi=200)
anim.save('field_sweep.mp4', dpi=200)

