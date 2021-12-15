import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from cycler import cycler
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl

mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14

# This function takes a 2D matrix T_mat and make a colormap of the values
def show_time_t(T_mat, title):
    plt.imshow(T_mat, cmap='viridis', extent=[0,1,0,1])
    plt.xlabel('x',fontsize=14)
    plt.ylabel('y',fontsize=14)
    cbar = plt.colorbar()
    cbar.set_label('p(x,y)', fontsize=14)
    plt.title(title, fontsize=16)
    plt.savefig(title+'.pdf', dpi=900)
    plt.show()

# This function taks two arrays y and x values and plot the results
def plot_prob(y, x, title):
    plt.plot(y, x)
    plt.xlabel("y", fontsize=14)
    plt.ylabel("Probability", fontsize=14)
    plt.title(title, fontsize=16)
    plt.grid()
    plt.tight_layout()
    plt.savefig(title+'.pdf', dpi=900)
    plt.show()



#Import files for problem 7 (one with wall and one without)
df1 = pd.read_csv("problem7_1.txt", sep=' ', header=None)
df2 = pd.read_csv("problem7_2.txt", sep=' ', header=None)


#Plot the deviation from 1
probabilities1 = np.array(df1.iloc[:,-1])
time1 = np.array(df1.iloc[:,-2])
smooth = savgol_filter(probabilities1, 51, 3)
plt.plot(time1, probabilities1)
plt.plot(time1, smooth)
plt.ylabel("Deviation from 1", fontsize=14)
plt.xlabel("Time", fontsize=14)
plt.title("Deviation without wall", fontsize=16)
plt.grid()
plt.savefig('Deviation_nowall.pdf', dpi=900)
plt.show()

#Plot the deviation from 1
probabilities2 = np.array(df2.iloc[:,-1])
time2 = np.array(df2.iloc[:,-2])
smooth = savgol_filter(probabilities2, 51, 3)
plt.plot(time2, probabilities2)
plt.plot(time2, smooth)
plt.ylabel("Deviation from 1", fontsize=14)
plt.xlabel("Time", fontsize=14)
plt.title("Deviation with wall and 2 slits", fontsize=16)
plt.grid()
plt.savefig('Deviation_wall.pdf', dpi=900)
plt.show()


#Problem 8
df8 = pd.read_csv("problem8.txt", sep=' ', header=None) #Import file
M = int(np.sqrt((len(df8.iloc[0,:])-2)/2)) #Calculate the M value


#T=0
T0 = np.array(df8.iloc[0,:-2])
T0imag = T0[1::2].reshape(M,M)
T0real = T0[::2].reshape(M,M)
print(T0imag.shape)
print(T0real.shape)
T0probs = T0real**2 + T0imag**2

show_time_t(T0imag, "Imag value at t=0")
show_time_t(T0real, "Real value at t=0")
show_time_t(T0probs, "2D probability at t=0")



#T=0.001
T1 = np.array(df8.iloc[40,:-2])
T1imag = T1[1::2].reshape((M,M))
T1real = T1[::2].reshape((M,M))
T1probs = T1imag**2 + T1real**2


show_time_t(T1imag, "Imag value at t=0.001")
show_time_t(T1real, "Real value at t=0.001")
show_time_t(T1probs, "2D probability at t=0.001")


#T=0.002
T2 = np.array(df8.iloc[-1,:-2])
T2imag = T2[1::2].reshape((M,M))
T2real = T2[::2].reshape((M,M))
T2probs = T2imag**2 + T2real**2


show_time_t(T2imag, "Imag value at t=0.002")
show_time_t(T2real, "Real value at t=0.002")
show_time_t(T2probs, "2D Probability at t=0.002")



#Problem 9

at08 = T2probs[:, int(0.8/0.005)]
at08 = at08/np.sum(at08)

y = np.linspace(0,1,201)
plot_prob(y, at08, "PDF of y when x = 0.8 with 2 slits")


#Problem 9.1

df91 = pd.read_csv("problem9_1slit.txt", sep=' ', header=None)
T91 = np.array(df91.iloc[-1,:-2])
T91imag = T91[1::2].reshape((M,M))
T91real = T91[::2].reshape((M,M))
T91probs = T91imag**2 + T91real**2


at08_91 = T91probs[:, int(0.8/0.005)]
at08_91 = at08_91/np.sum(at08_91)

y = np.linspace(0,1,201)
plot_prob(y, at08_91, "PDF of y when x = 0.8 with 1 slit")


#Problem 9.3

df92 = pd.read_csv("problem9_3slit.txt", sep=' ', header=None)
T92 = np.array(df92.iloc[-1,:-2])
T92imag = T92[1::2].reshape((M,M))
T92real = T92[::2].reshape((M,M))
T92probs = T92imag**2 + T92real**2


at08_92 = T92probs[:, int(0.8/0.005)]
at08_92 = at08_92/np.sum(at08_92)

y = np.linspace(0,1,201)
plot_prob(y, at08_92, "PDF of y when x = 0.8 with 3 slits")
