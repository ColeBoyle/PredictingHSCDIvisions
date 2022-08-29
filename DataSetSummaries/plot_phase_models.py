import matplotlib.pyplot as plt
import data
import matplotlib.transforms as mtransforms
import fit_pwlf as fit

# Plot 1 to 6 phase models on data in single figure

# Set data to use
x = data.AubertData_age
y = data.AubertData_gran

# Do all n-phase model fits (1 to 6 phases)
fits = [fit.fit_pwlf(1, 1, x, y),
        fit.fit_pwlf(2, 1, x, y),
        fit.fit_pwlf(3, 1, x, y),
        fit.fit_pwlf(4, 1, x, y),
        fit.fit_pwlf(5, 1, x, y),
        fit.fit_pwlf(6, 1, x, y)]

# Plot the fit of each phase model on data
fontsize= 7
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
fig, axs = plt.subplot_mosaic([['A', 'B'],
                               ["C", "D"],
                               ["E", "F"]],
                              figsize=(155*mm, 120*mm), constrained_layout=True)

for label, ax in axs.items():
    trans = mtransforms.ScaledTranslation(10 / 72, -5 / 72, fig.dpi_scale_trans)
    ax.text(-0.05, 1.2, label, transform=ax.transAxes + trans,
            fontweight="bold", fontsize=fontsize+2, verticalalignment='top')

for i, ax in enumerate(axs.values()):
    ax.plot(x, y, ".", color="grey")
    ax.plot(fits[i].xHat, fits[i].yHat, color="blue")
    ax.legend(title=f"    {i + 1}-phase\n$R^2 = {fits[i].pwlf.r_squared():.4f}$\n"
                    f"$\Delta$ AIC = {fits[i].aic - fits[-1].aic:.2f}", shadow=True)

axs.get("E").set_xlabel("Age (years)")
axs.get("F").set_xlabel("Age (years)")
fig.supylabel("Telomere Length (kbp)")

plt.savefig('models.pdf', dpi=300, format="PDF", pad_inches=0)
plt.show()
