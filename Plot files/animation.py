import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
from matplotlib.animation import FuncAnimation
import pandas as pd

df8 = pd.read_csv("problem8.txt", sep=' ', header=None)
M = int(np.sqrt((len(df8.iloc[0,:])-2)/2))

z = []
for i in range(len(df8)):
    T = np.array(df8.iloc[i,:-2])
    Timag = T[1::2].reshape(M,M)
    Treal = T[::2].reshape(M,M)
    Tprobs = Treal**2 + Timag**2
    z.append(Tprobs)

xmin, ymin = 0,0
xmax, ymax = 1,1
t0 = 0

dt = 2.5e-5

def animate(z, name):
    z = np.sqrt(z) # This is added just to make the probabilities even more visibile in the amination.
    fig = plt.figure()
    ax = plt.gca()

    norm = mat.cm.colors.Normalize(vmin=0.0, vmax=np.max(z[0]))
    img = ax.imshow(z[0], extent=[xmin, xmax, ymin, ymax], cmap=plt.get_cmap('viridis'), norm=norm)

    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label('P(x,y)', fontsize=12)
    cbar.ax.tick_params(labelsize=12)


    time_txt = plt.text(0.95,0.95, 't={:.3e}'.format(t0), color='white', horizontalalignment='right', verticalalignment='top', fontsize=12)

    def animation(i):
        norm = mat.cm.colors.Normalize(vmin=0.0, vmax=np.max(z[i]))
        img.set_norm(norm)

        img.set_data(z[i])

        current_time = t0 + i * dt
        time_txt.set_text('t={:.3e}'.format(current_time))

        return img


    anim = FuncAnimation(fig,animation, interval=1, frames=np.arange(0,len(z), 2), repeat=False, blit=0)

    plt.show()

    anim.save('./animation_'+name+'.gif',dpi=700, writer='ffmpeg', bitrate=-1, fps=10)

animate(z, '8')
