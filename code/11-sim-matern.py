

import numpy as np
import matplotlib.pyplot as plt; import warnings; warnings.filterwarnings("ignore")
plt.rcParams.update(**{'font.size': 12, 'axes.labelsize': 12, 'legend.fontsize': 12,
'figure.dpi': 400, 'axes.titlesize': 12, 'xtick.labelsize': 12, 'ytick.labelsize': 12})
import seaborn as sns; 
sns.set_style("whitegrid", {"grid.color": ".9"})
sns.set_palette("tab10")
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern



# RdPu, but I change the end chroma to 100
YlRdPu = ["#490062", "#963E87", "#E775A8", "#FFAEAD", "#FFE08C"]
from matplotlib.colors import LinearSegmentedColormap
clist = sns.color_palette(YlRdPu)
my_cmap = LinearSegmentedColormap.from_list('my_cmap',clist,80)
sns.set_palette(YlRdPu)



### PLOT 1D SIMULATIONS OF GAUSSIAN PROCESS WITH MATERN COVARIANCE
nu_array = np.array([1.0,2.0,4.0])
n_samples = 5
x = np.linspace(0, 6, 2001)
X = x.reshape(-1, 1)



fig, axs = plt.subplots(ncols=np.size(nu_array), sharex=False, sharey=False, figsize=(12, 4))
for (i,nu) in enumerate(nu_array):
    
    kernel = 1.0 * Matern(length_scale=1.0, length_scale_bounds=(1e-1, 10.0), nu=nu)
    gpr = GaussianProcessRegressor(kernel=kernel, random_state=100)
    y_samples = gpr.sample_y(X, n_samples)
    for n in np.arange(n_samples):
        sns.lineplot(x=x,y=y_samples[:,n],ax=axs[i],linewidth=1,alpha=0.9)
    axs[i].set_xlim([0,6])
    axs[i].set_ylim([-3,3])
    axs[i].set_xticks([0,1,2,3,4,5,6])
    axs[i].set_aspect("equal") 
    
    
    
    
    
    
### PLOT 2D SIMULATIONS OF GAUSSIAN PROCESS WITH MATERN COVARIANCE
from matplotlib import cm
from matplotlib.ticker import LinearLocator
plt.rcParams.update(**{'font.size': 11, 'axes.labelsize': 11, 'legend.fontsize': 11,
'figure.dpi': 400, 'axes.titlesize': 11, 'xtick.labelsize': 11, 'ytick.labelsize': 11})

nu_array = np.array([1.0,2.0,4.0])
x = np.linspace(0, 5, 43)
y = np.linspace(0, 5, 43)
xv, yv = np.meshgrid(x, y)
X = np.concatenate( ( xv.flatten().reshape(-1,1) , yv.flatten().reshape(-1,1) ), axis=1 )
    
fig, axs = plt.subplots(ncols=np.size(nu_array), subplot_kw={"projection": "3d"}, sharex=False, sharey=False, figsize=(12, 4))
for (i,nu) in enumerate(nu_array):

    ax = axs[i]
    kernel = 1.0 * Matern(length_scale=1.0, length_scale_bounds=(1e-1, 10.0), nu=nu)
    gpr = GaussianProcessRegressor(kernel=kernel, random_state=0)
    z_samples = gpr.sample_y(X, 2)[:,1].reshape(-1,1)
    zv = np.zeros_like(xv)
    nx = np.size(x)
    for n in np.arange(nx):
        zv[n,:] = np.squeeze(z_samples[n*nx:(n+1)*nx])
    
    surf = ax.plot_surface(xv, yv, zv, cmap=my_cmap,linewidth=0, antialiased=False)
    ax.zaxis.set_major_locator(LinearLocator(6))
    ax.zaxis.set_major_formatter('{x:.0f}')
    ax.set_xlim([0,5])
    ax.set_ylim([0,5])
    ax.set_zlim([-3,4])
    ax.set_zticks([-2,0,2,4])
    fig.colorbar(surf, shrink=0.5, aspect=10, pad=0.1)


