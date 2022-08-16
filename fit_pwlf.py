import numpy as np
import pwlf

# Use pwlf to fit an optimal n-interval piecewise linear model to telomere length data, then recover division estimates
#  from fitted model.

class fit_pwlf:

    def __init__(self, interval, degree, xdata, ydata):
        self.n_intervals = interval
        self.degree = degree
        self.statStr = ""
        self.xdata = xdata
        self.ydata = ydata
        self.aic = 0
        self.pwlf = pwlf.PiecewiseLinFit(self.xdata, self.ydata, seed=111, degree=self.degree)
        breakpoints = self.pwlf.fit(self.n_intervals)
        self.r_values_50bp = []
        self.r_values_150bp = []
        self.r_values_100bp = []
        self.r_values_30bp = []
        self.xHat = np.linspace(min(self.xdata), max(self.xdata), num=10000)
        self.yHat = self.pwlf.predict(self.xHat)

        slopes = self.pwlf.slopes
        self.aic = aic(self.pwlf.n_parameters + 1, len(self.xdata), self.pwlf.ssr)

        self.statStr = f"$R^2$ = {self.pwlf.r_squared():.4f}" \
                       f"\nAIC = {self.aic:.2f}"

        self.r_values_50bp = slopes / (-0.05)
        self.r_values_150bp = slopes / (-0.15)
        self.r_values_100bp = slopes / (-0.1)
        self.r_values_30bp = slopes / (-0.03)

        # Determine r(t):
        rt_50bp = []
        rt_150bp = []
        rt_100bp = []
        rt_30bp = []
        time = []
        for i in range(len(breakpoints) - 1):
            tp = np.linspace(breakpoints[i] + 0.0001, breakpoints[i + 1])
            rt_100bp += (self.r_values_100bp[i] * np.ones(len(tp))).tolist()
            rt_150bp += (self.r_values_150bp[i] * np.ones(len(tp))).tolist()
            rt_50bp += (self.r_values_50bp[i] * np.ones(len(tp))).tolist()
            rt_30bp += (self.r_values_30bp[i] * np.ones(len(tp))).tolist()
            time += tp.tolist()

        self.rt_100bp = np.array(rt_100bp)
        self.rt_50bp = np.array(rt_50bp)
        self.rt_150bp = np.array(rt_150bp)
        self.rt_30bp = np.array((rt_30bp))
        self.time = np.array(time)

        def r(t, dL):
            if dL == 50:
                R = self.r_values_50bp
            elif dL == 100:
                R = self.r_values_100bp
            elif dL == 150:
                R = self.r_values_150bp
            else:
                R = self.r_values_30bp

            return R

        def Dbar(t, dL):
            result = [0]
            q = 0
            for i, _ in enumerate(breakpoints[1:]):
                q += r(t, dL)[i] * (_ - breakpoints[i])
                result += [q]
            return np.array(result)

        self.res = breakpoints
        self.div30 = Dbar(time, 30)
        self.div50 = Dbar(time, 50)
        self.div100 = Dbar(time, 100)

        # Print Results:
        print("-" * 16, f"{self.pwlf.n_parameters - 1}-phase model: Data Set:{len(self.xdata)}", "-" * 16)
        print("|  Age Range     | Division Rate Range | Number of Divisions |")
        for i, _ in enumerate(breakpoints[1:]):
            print("-" * 62)
            print(
                f"| {breakpoints[i]:.2f} to {_:.2f} | {self.r_values_100bp[i]:.4f} to {self.r_values_30bp[i]:.4f}       | {Dbar(time, 100)[i + 1] - Dbar(time, 100)[i]:.2f} to {Dbar(time, 30)[i + 1] - Dbar(time, 30)[i]:.2f}        |")
        print("-" * 62)
        print(
            f"| Total Divisions: {Dbar(time, 100)[-1]:.2f} to {Dbar(time, 30)[-1]:.2f} Over {np.max(self.xdata) - np.min(self.xdata)} years           |")
        print("-" * 62)


def aic(k, n, sse):
    sigma_hat_sq = sse / n
    LL = - n / 2 * np.log(sigma_hat_sq) - 1 / (2 * sigma_hat_sq) * sse - n / 2 * np.log(2 * np.pi)
    return 2 * k - 2 * LL
