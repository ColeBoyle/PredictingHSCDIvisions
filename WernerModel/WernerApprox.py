import numpy as np
import matplotlib.pyplot as plt
import data


dcr_N0 = 0.068
c = 10.2
p = 0.44

n = 35
N0 = 100

dc = 0.065
r = (dcr_N0 * N0)/dc
T = 102
t = np.linspace(0, T, 3000)


plt.figure()
plt.plot(data.AubertData_age, data.AubertData_gran, ".", color="grey")
plt.plot(t, c - dc*(1+p)/p * np.log(p*dcr_N0/dc *t + 1), label=r"$\bar L(t) = c - \Delta c \frac{1+p}{p}\log(p\frac{r}{N_0} t + 1)$")
plt.plot(t, c - (1+p)/p * np.log(p*dcr_N0 *t + 1), label = r"$\bar L'(t) = c - \frac{1+p}{p}\log(p\frac{r \Delta c }{N_0} t + 1)$")
plt.xlabel("Time (t) (years)")
plt.ylabel("Telomere Length (kbp)")
plt.legend()
plt.savefig("wernerApprox.pdf", dpi=300, format="PDF", pad_inches=0)

plt.show()