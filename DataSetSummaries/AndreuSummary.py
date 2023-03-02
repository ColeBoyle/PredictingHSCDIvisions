import fit_pwlf as fit
import numpy as np
import data
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec

# Make figure summarizing results for Andreu et al. data

x = data.AndreuData_age
y = data.AndreuData_gran


twoPhase = fit.fit_pwlf(2, 1, x, y, [18, 34, 82])
threePhase = fit.fit_pwlf(1, 1, x, y, [18, 82])

# Set up figure
fontsize = 7
plt.rcParams.update({'font.size': fontsize})
plt.rcParams.update({'figure.dpi': 300})
plt.rcParams.update({'savefig.dpi': 300})
plt.rcParams.update({'lines.markersize': 2.5})
plt.rcParams.update({'lines.linewidth': 1.0})
plt.rcParams.update({'ytick.major.width': 0.5})
plt.rcParams.update({'xtick.major.width': 0.5})
plt.rcParams.update({'ytick.major.size': 2})
plt.rcParams.update({'xtick.major.size': 2})
plt.rcParams.update({'ytick.minor.size': 1.5})
plt.rcParams.update({'xtick.minor.size': 1.5})
plt.rcParams.update({'legend.edgecolor': 'black'})

mm = 0.03937
fig = plt.figure(figsize=(155 * 1 / 2 * mm, 115 * mm))
gs = gridspec.GridSpec(3, 1)
ax0 = fig.add_subplot(gs[0])
#ax1 = fig.add_subplot(gs[1], sharey=ax0)

ax4 = fig.add_subplot(gs[1], sharex=ax0)
#ax5 = fig.add_subplot(gs[3], sharey=ax4, sharex=ax1)

ax2 = fig.add_subplot(gs[2])
#ax3 = fig.add_subplot(gs[5], sharey=ax2)

#plt.setp(ax1.get_yticklabels(), visible=False)
#plt.setp(ax3.get_yticklabels(), visible=False)
#plt.setp(ax5.get_yticklabels(), visible=False)

plt.setp(ax0.get_xticklabels(), visible=False)
#plt.setp(ax1.get_xticklabels(), visible=False)

ax4.set_xticks([18, 40, 60, 80, 100])
#ax5.set_xticks([23, 40, 60, 80, 100])

# Plot model on Data
ax0.text(-0.15, 1.1, "A1", transform=ax0.transAxes,
         size=fontsize + 2, weight='bold')
ax0.set_title("1-phase Model", weight='bold')
ax0.set_ylabel("Telomere Length (kbp)")
ax0.plot(x, y, ".", color="grey", markersize=1)

ax0.set_ylim(3.5, 15)
ax0.plot(threePhase.xHat, threePhase.yHat, "blue")
ax0.legend(title=f"$R^2= {threePhase.pwlf.r_squared():.4f}$")
#                f"$\Delta$ AIC = 0")

#ax1.set_title("2-phase Model", weight='bold')
#ax1.plot(x, y, ".", color="grey", markersize=1)
#ax1.plot(twoPhase.xHat, twoPhase.yHat, "blue")
#ax1.plot([23, 23], [0, threePhase.pwlf.predict(23)], "k--")
#ax1.plot(23, threePhase.pwlf.predict(23), "ko")
#ax1.legend(title=f"$R^2= {twoPhase.pwlf.r_squared():.4f}$\n"
#                 f"$\Delta$ AIC = {twoPhase.aic - threePhase.aic:.2f}")

# Division Rates
ax2.text(-0.15, 1.05, "C1", transform=ax2.transAxes,
         size=fontsize + 2, weight='bold')
ax2.set_ylabel("Division Rate (1/year)")

phases_2 = np.array(['18-81 yrs'])
x_axis_2 = np.array([1])
ax2.bar(x_axis_2 + 0.1, threePhase.r_values_100bp, width=0.1, color="#2c7bb6", label="100 bp", zorder=3)
ax2.bar(x_axis_2, threePhase.r_values_65bp, width=0.1, color="black", label="65 bp", zorder=3)
ax2.bar(x_axis_2 - 0.1, threePhase.r_values_30bp, width=0.1, color="#d7191c", label="30 bp", zorder=3)
ax2.set_xticks(ticks=x_axis_2, labels=phases_2)
ax2.tick_params(which="both")
#ax3.tick_params(which="both")
#ax3.grid(which="both", zorder=0)
ax2.grid(which="both", zorder=0)
ax2.set_xlabel("Phases")

phases_3 = np.array(['18-23 yrs', '23-81 yrs'])
x_axis_3 = np.array([1, 1.35])
#ax3.bar(x_axis_3 + 0.1, twoPhase.r_values_100bp, width=0.1, color="#2c7bb6", label="100 bp", zorder=3)
#ax3.bar(x_axis_3, twoPhase.r_values_65bp, width=0.1, color="black", label="65 bp", zorder=3)
#ax3.bar(x_axis_3 - 0.1, twoPhase.r_values_30bp, width=0.1, color="#d7191c", label="30 bp", zorder=3)
#ax3.set_xticks(ticks=x_axis_3, labels=phases_3)
#ax3.set_xlabel("Phases")

ax2.set_yscale('log', subs=[2, 3, 4, 5, 6, 7, 8, 9])
ax2.legend(title=f"$\Delta L$")
ax2.set_yticks([1])
ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

# Division Numbers
ax4.text(-0.15, 1.1, "B1", transform=ax4.transAxes,
         size=fontsize + 2, weight='bold')
ax4.set_ylabel("Divisions")
ax4.set_xlabel("Age (years)")
ax4.grid()
ax4.plot(threePhase.breakpoints, threePhase.div_30bp, color="#d7191c")
ax4.plot(threePhase.breakpoints, threePhase.div_65bp, color="black", linewidth=1.5)
ax4.plot(threePhase.breakpoints, threePhase.div_100bp, color="#2c7bb6")
ax4.set_ylim(0, 62)

#ax5.plot([23, 23], [0, 160], "k--")
#
#ax5.grid()
#ax5.plot(twoPhase.breakpoints, twoPhase.div_30bp, color="#d7191c")
#ax5.plot(twoPhase.breakpoints, twoPhase.div_65bp, color="black", linewidth=1.5)
#ax5.plot(twoPhase.breakpoints, twoPhase.div_100bp, color="#2c7bb6")
#ax5.set_xlabel("Age (years)")

gs.tight_layout(fig, h_pad=0.2, w_pad=-0.4)
plt.savefig('AndreuSummary.pdf', dpi=300, format="PDF", pad_inches=0)

plt.show()
