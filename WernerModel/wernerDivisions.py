import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Obtain division estimates from Werner et al. (2015) model

# Reported parameters
c = 10.2
p = 0.44
rdc_N0 = 0.068

# Additional parameters
N0 = 100
n = 35
T = 102  # Maximum time

# list of dc values to use
dc_list = [0.1, 0.065, 0.03]

# list of division results
Dbar_list = []

# Initial condition
x0 = [N0] + [0 for i in range(n)]

# Times to sample
t = np.linspace(0, T, 300)

def wernerODE(t, z, r):
    N = N0 + r * p * t
    return [-r / N * z[0]] + [r / N * ((1 + p) * z[i - 1] - z[i])
                              for i in range(1, n)] + [r / N * z[n - 1]]

# Solve ODE system for the different dc values
for dc in dc_list:
    r = rdc_N0 * N0 / dc
    sol = solve_ivp(wernerODE, [0, T], x0, args=(r,), dense_output=True)
    A = sol.sol(t).T

    Dbar_list += [np.sum(np.matmul(A, np.diag(np.arange(0, n + 1, 1))), axis=1) / np.sum(A, axis=1)]

# Plot results
fontsize = 7
plt.rcParams.update({'font.size': fontsize})
plt.rcParams.update({'figure.dpi': 300})
plt.rcParams.update({'savefig.dpi': 300})
plt.rcParams.update({'lines.markersize': 1})
plt.rcParams.update({'lines.linewidth': 1.5})
plt.rcParams.update({'ytick.major.width': 0.5})
plt.rcParams.update({'xtick.major.width': 0.5})
plt.rcParams.update({'figure.subplot.bottom': 0.15})
plt.rcParams.update({'figure.subplot.left': 0.15})
plt.rcParams.update({'ytick.major.size': 2})
plt.rcParams.update({'xtick.major.size': 2})
plt.rcParams.update({'ytick.minor.size': 1.5})
plt.rcParams.update({'xtick.minor.size': 1.5})
plt.rcParams.update({'legend.edgecolor': 'black'})
mm = 0.03937

fig, ax1 = plt.subplots(figsize=(155*1/2*mm, 60*mm))
ax1.plot(t, Dbar_list[0], label="100 bp", color="#d7191c")
ax1.plot(t, Dbar_list[1], label="65 bp", color="black", linewidth=2)
ax1.plot(t, Dbar_list[2], label="30 bp", color="#2c7bb6")
ax1.legend(title=r"$\Delta L$")
ax1.set_ylabel("Divisions")
ax1.set_xlabel("Age (years)")
ax1.grid()
plt.savefig('wernerDivision.pdf', dpi=300, format="PDF", pad_inches=0)
plt.show()

# Repeat with predicted telomere length function

def logDbar(t, params):
    p, r_N0 = params
    return (1+p)/p*np.log(p*r_N0*t + 1)

logDbar_list = []

for dc in dc_list:
    logDbar_list += [logDbar(t, [p, rdc_N0/dc])]

fig2, ax2 = plt.subplots(figsize=(155*1/2*mm, 80*mm ))
ax2.plot(t, logDbar_list[0], label="100 bp", color="#d7191c")
ax2.plot(t, logDbar_list[1], label="65 bp", color="black", linewidth=2)
ax2.plot(t, logDbar_list[2], label="30 bp", color="#2c7bb6")
ax2.legend(title=r"$\Delta L$")
ax2.set_ylabel("Divisions")
ax2.set_title(r"Predicted divisons from logarithmic $\overline{D}(t)$")
ax2.set_xlabel("Age (years)")
ax2.grid()
#plt.show()
