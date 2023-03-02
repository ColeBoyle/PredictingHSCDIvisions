import fit_pwlf as fit
import data
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec
import string
import numpy as np

# Make figure summarizing results for Aubert et al. data

x = data.AubertData_age
y = data.AubertData_gran


bp_2 = [0, 20.9, 102]
bp_3 = [0, 0.9510914, 24.2078681, 102]

div_2_30bp = [0.27091, 0.08842]
div_2_65bp = [1.2503538, 0.4080923]
div_2_100bp = [0.81273, 0.26526]


div_3_30bp = [19.0286667, 2.0854667, 0.8844667]
div_3_65bp = [8.7824615, 0.9625231, 0.4082154]
div_3_100bp = [5.70860, 0.62564, 0.26534]

twoPhase = fit.fit_pwlf(2, 1, x, y, bp_2)
threePhase = fit.fit_pwlf(3, 1, x, y, bp_3)

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
fig = plt.figure(figsize=(155 * mm, 115 * mm))
gs = gridspec.GridSpec(3, 2)

ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1], sharey=ax0)

ax2 = fig.add_subplot(gs[4])
ax3 = fig.add_subplot(gs[5], sharey=ax2)

ax4 = fig.add_subplot(gs[2], sharex=ax0)
ax5 = fig.add_subplot(gs[3], sharey=ax4, sharex=ax1)

plt.setp(ax1.get_yticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)

plt.setp(ax0.get_xticklabels(), visible=False)
plt.setp(ax1.get_xticklabels(), visible=False)

ax4.set_xticks([0.95, 24.2, 40, 60, 80, 100])
ax5.set_xticks([0, 20.9, 40, 60, 80, 100])

# Fit models to data
ax0.text(-0.15, 1.1, string.ascii_uppercase[0], transform=ax0.transAxes,
         size=fontsize + 2, weight='bold')
ax0.set_title("3-phase Model", weight='bold')
ax0.set_ylabel("Telomere Length (kbp)")
ax0.plot(x, y, ".", color="grey", markersize=1)
ax0.plot([bp_3[2], bp_3[2]], [0, threePhase.pwlf.predict(bp_3[2])], "k--")
ax0.plot(bp_3[2], threePhase.pwlf.predict(bp_3[2]), "ko")
ax0.plot(bp_3[1], threePhase.pwlf.predict(bp_3[1]), "ko")
ax0.plot([bp_3[1], bp_3[1]], [0, threePhase.pwlf.predict(bp_3[1])], "k--")
ax0.set_ylim(3.5, 15)
ax0.plot(threePhase.xHat, threePhase.yHat, "blue")
ax0.legend(title=f"$R^2= {threePhase.pwlf.r_squared():.4f}$\n"
                 f"$\Delta$ AIC = 0", shadow=True)

ax1.set_title("2-phase Model", weight='bold')
ax1.plot(x, y, ".", color="grey", markersize=1)
ax1.plot(twoPhase.xHat, twoPhase.yHat, "blue")
ax1.plot([bp_2[1], bp_2[1]], [0, threePhase.pwlf.predict(bp_2[1])], "k--")  # change
ax1.plot(bp_2[1], threePhase.pwlf.predict(bp_2[1]), "ko")  # change
ax1.legend(title=f"$R^2= {twoPhase.pwlf.r_squared():.4f}$\n"
                 f"$\Delta$ AIC = {twoPhase.aic - threePhase.aic:.2f}", shadow=True)

# Division Rates
ax2.text(-0.15, 1.1, string.ascii_uppercase[2], transform=ax2.transAxes,
         size=fontsize + 2, weight='bold')
ax2.set_ylabel("Division Rate (1/year)")

phases_3 = np.array(['0-0.95 yrs', '0.95-24.2 yrs', '24.2-102 yrs'])
x_axis_3 = np.array([1, 1.5, 2])
ax2.bar(x_axis_3 + 0.1, threePhase.r_values_100bp, width=0.1, color="#2c7bb6", label="100 bp", zorder=3)
ax2.bar(x_axis_3, threePhase.r_values_65bp, width=0.1, color="black", label="65 bp", zorder=3)
ax2.bar(x_axis_3 - 0.1, threePhase.r_values_30bp, width=0.1, color="#d7191c", label="30 bp", zorder=3)
ax2.set_xticks(ticks=x_axis_3, labels=phases_3)
ax2.tick_params(which="both")
ax3.tick_params(which="both")
ax3.grid(which="both", zorder=0)
ax2.grid(which="both", zorder=0)
ax2.set_xlabel("Phases")
ax2.set_yscale('log')
ax2.set_yticks([1, 10])
ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

phases_3 = np.array([f'0-{bp_2[1]:.1f} yrs', f'{bp_2[1]:.1f}-102 yrs'])
x_axis_3 = np.array([1, 1.4])
ax3.bar(x_axis_3 + 0.1, twoPhase.r_values_100bp, width=0.1, color="#2c7bb6", label="100 bp", zorder=3)
ax3.bar(x_axis_3, twoPhase.r_values_65bp, width=0.1, color="black", label="65 bp", zorder=3)
ax3.bar(x_axis_3 - 0.1, twoPhase.r_values_30bp, width=0.1, color="#d7191c", label="30 bp", zorder=3)
ax3.set_xticks(ticks=x_axis_3, labels=phases_3)
ax3.legend(title=f"$\Delta L$", shadow=True)
ax3.set_xlabel("Phases")

# Division Numbers
ax4.text(-0.15, 1.1, string.ascii_uppercase[1], transform=ax4.transAxes,
         size=fontsize + 2, weight='bold')
ax4.set_ylabel("Divisions")
ax4.set_xlabel("Age (years)")
ax4.grid()
ax4.plot(threePhase.breakpoints, threePhase.div_30bp, color="#d7191c")
ax4.plot(threePhase.breakpoints, threePhase.div_65bp, color="black", linewidth=1.5)
ax4.plot(threePhase.breakpoints, threePhase.div_100bp, color="#2c7bb6")
ax4.plot([bp_3[2], bp_3[2]], [0, 160], "k--")
ax4.plot([bp_3[1], bp_3[1]], [0, 160], "k--")
ax4.set_ylim(0, 153)

ax5.plot([bp_2[1], bp_2[1]], [0, 160], "k--")
ax5.grid()
ax5.plot(twoPhase.breakpoints, twoPhase.div_30bp, color="#d7191c")
ax5.plot(twoPhase.breakpoints, twoPhase.div_65bp, color="black", linewidth=1.5)
ax5.plot(twoPhase.breakpoints, twoPhase.div_100bp, color="#2c7bb6")
ax5.set_xlabel("Age (years)")

gs.tight_layout(fig, h_pad=0, w_pad=-0.5)
plt.savefig('AubertSummary.pdf', dpi=300, format="PDF", pad_inches=0)
plt.show()