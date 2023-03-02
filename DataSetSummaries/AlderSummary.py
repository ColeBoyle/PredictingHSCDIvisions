import fit_pwlf as fit
import data
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import matplotlib

# Make figure summarizing results for Alder et al. data

x = data.AlderData_age
y = data.AlderData_gran

twoPhase = fit.fit_pwlf(2, 1, x, y, [0, 26, 82])

# Set up figure
fontsize = 7
plt.rcParams.update({'font.size': fontsize}) # textwidth 156 mm
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
fig = plt.figure(figsize=(155*1/2*mm, 115*mm))
gs = gridspec.GridSpec(3, 1)

ax0 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[2])
ax4 = fig.add_subplot(gs[1], sharex=ax0)

plt.setp(ax0.get_xticklabels(), visible=False)

ax4.set_xticks([0, 26, 40, 60, 80, 100])

# Fit model to data
ax0.text(-0.15, 1.1, "A2", transform=ax0.transAxes,
         size=fontsize+2, weight='bold')
ax0.set_title("2-phase Model", weight='bold')
ax0.set_ylabel("Telomere Length (kbp)")
ax0.plot(x, y, ".", color="grey", markersize=1)
ax0.plot([26,26], [0, twoPhase.pwlf.predict(26)], "k--")
ax0.plot(26, twoPhase.pwlf.predict(26), "ko")
ax0.set_ylim(3.5, 12)
ax0.plot(twoPhase.xHat, twoPhase.yHat, "blue", markersize=1)

ax0.legend(title=f"$R^2= {twoPhase.pwlf.r_squared():.4f}$")

# Division Rates
ax2.text(-0.15, 1.1, "C2", transform=ax2.transAxes,
         size=fontsize+2, weight='bold')
ax2.set_ylabel("Division Rate (1/year)")


phases_3 = np.array(['0-26 yrs', '26-82 yrs'])
x_axis_3 = np.array([1, 1.35])
ax2.bar(x_axis_3 + 0.1, twoPhase.r_values_100bp, width=0.1, color="#2c7bb6", label="100 bp", zorder=3)
ax2.bar(x_axis_3, twoPhase.r_values_65bp, width=0.1, color="black", label="65 bp", zorder=3)
ax2.bar(x_axis_3 - 0.1, twoPhase.r_values_30bp, width=0.1, color="#d7191c", label="30 bp", zorder=3)
ax2.set_xticks(ticks=x_axis_3, labels=phases_3)
ax2.tick_params(which="both")

ax2.grid(which="both", zorder=0)
ax2.set_xlabel("Phases")
ax2.set_yscale('log')
ax2.set_yticks([1])
ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

# Division Numbers
ax4.text(-0.15, 1.1, "B2", transform=ax4.transAxes,
         size=fontsize+2, weight='bold')
ax4.set_ylabel("Divisions")
ax4.set_xlabel("Age (years)")
ax4.grid()
ax4.plot(twoPhase.breakpoints, twoPhase.div_30bp, color="#d7191c")
ax4.plot(twoPhase.breakpoints, twoPhase.div_65bp, color="black", linewidth=1.5)
ax4.plot(twoPhase.breakpoints, twoPhase.div_100bp, color="#2c7bb6")
ax4.plot([26, 26],[0, 160],"k--")
ax4.set_ylim(0, 120)

gs.tight_layout(fig, h_pad=0, w_pad=-0.5)
plt.savefig('AlderSummary.pdf', dpi=300, format="PDF", pad_inches=0)
plt.show()
