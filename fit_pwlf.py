import numpy as np
import pwlf


# Use pwlf to fit an optimal n-interval piecewise linear model to telomere length data, then recover division estimates
#  from fitted model.

class fit_pwlf:

    def __init__(self, interval, degree, xdata, ydata, bp):
        self.bp = bp
        self.n_intervals = interval
        self.degree = degree
        self.statStr = ""
        self.xdata = xdata
        self.ydata = ydata
        self.aic = 0
        self.pwlf = pwlf.PiecewiseLinFit(self.xdata, self.ydata, seed=111, degree=self.degree)
        #print(self.pwlf.fit(self.n_intervals))
        self.breakpoints = bp
        self.pwlf.fit_with_breaks(self.bp)
        self.xHat = np.linspace(min(self.xdata), max(self.xdata), num=10000)
        self.yHat = self.pwlf.predict(self.xHat)
        self.slopes = self.pwlf.slopes
        self.aic = aic(self.pwlf.n_parameters + 1 + len(bp), len(self.xdata), self.pwlf.ssr)
        self.bic = bic(self.pwlf.n_parameters + 1, len(self.xdata), self.pwlf.ssr)
        self.statStr = f"$R^2$ = {self.pwlf.r_squared():.4f}" \
                       f"\nAIC = {self.aic:.2f}"
        self.r_values_100bp, self.rt_100bp, self.div_100bp = self.getDiv(0.1)
        self.r_values_65bp, self.rt_50bp, self.div_50bp = self.getDiv(0.05)
        self.r_values_150bp, self.rt_150bp, self.div_150bp = self.getDiv(0.15)
        self.r_values_30bp, self.rt_30bp, self.div_30bp = self.getDiv(0.03)
        self.r_values_65bp, self.rt_65bp, self.div_65bp = self.getDiv(0.065)

        # Print Results:
        print("-" * 16, f"{self.pwlf.n_parameters - 1}-phase model: Data Set:{len(self.xdata)}", "-" * 16)
        print("|  Age Range     | Division Rate Range | Number of Divisions |")
        for i, _ in enumerate(self.breakpoints[1:]):
            print("-" * 62)
            print(
                f"| {self.breakpoints[i]:.2f} to {_:.2f} | {self.r_values_65bp[i]:.4f} ({self.r_values_100bp[i]:.4f} to {self.r_values_30bp[i]:.4f})"
                f"       | {self.div_65bp[i + 1] - self.div_65bp[i]:.4f} ({self.div_100bp[i + 1] - self.div_100bp[i]:.4f} to {self.div_30bp[i + 1] - self.div_30bp[i]:.4f})        "
                f"|")
        print("-" * 62)
        print(f"| Total Divisions: {self.div_65bp[-1]:.4f} ({self.div_100bp[-1]:.4f} to {self.div_30bp[-1]:.4f}) Over "
              f"{np.max(self.xdata) - np.min(self.xdata)} years           |")
        print("-" * 62)

    def getDiv(self, bp_per_div):
        r_vals = self.slopes / (-bp_per_div)

        # Determine r(t):
        rt = []
        time = []
        for i in range(len(self.breakpoints) - 1):
            tp = np.linspace(self.breakpoints[i] + 0.0001, self.breakpoints[i + 1])
            rt += (r_vals[i] * np.ones(len(tp))).tolist()
            time += tp.tolist()
        rt = np.array(rt)
        self.time = np.array(time)

        result = [0]
        q = 0
        for i, _ in enumerate(self.breakpoints[1:]):
            q += r_vals[i] * (_ - self.breakpoints[i])
            result += [q]
        Dbar = np.array(result)

        return [r_vals, rt, Dbar]


def aic(k, n, sse):
    sigma_hat_sq = sse / n
    LL = - n / 2 * np.log(sigma_hat_sq) - 1 / (2 * sigma_hat_sq) * sse - n / 2 * np.log(2 * np.pi)
    result = 2 * k - 2 * LL
    return result


def bic(k, n, sse):
    sigma_hat_sq = sse / n
    LL = - n / 2 * np.log(sigma_hat_sq) - 1 / (2 * sigma_hat_sq) * sse - n / 2 * np.log(2 * np.pi)
    result = np.log(n) * k - 2 * LL
    return result
