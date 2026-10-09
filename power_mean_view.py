"""Tần số dao động tự do như một trung bình lũy thừa của độ cứng cát tuyến năng lượng.

Cho  x'' + g(x) = 0,  g lẻ, x(0) = A,  U(x) = int_0^x g. Đổi biến x = A sin(theta), theta đều trên [0, pi/2]:

    K(theta) = 2 [U(A) - U(A sin theta)] / (A^2 cos^2 theta)      (độ cứng cát tuyến năng lượng)
    1/omega  = E[K^(-1/2)],   tức  omega^2 = M_{-1/2}(K)            (chính xác)
    E[K]     = (2/A) <g(A cos t) cos t>  = M_1(K)                   (tuyến tính hóa cổ điển = cân bằng điều hòa bậc 1)

Bất đẳng thức trung bình lũy thừa cho, với mọi g:
    M_{-1}(K) <= omega^2 <= M_0(K) <= M_1(K).
Khai triển theo độ phân tán của K cho ln M_r ~ tuyến tính theo r, nên
    omega^2 ~ M_1^(1/4) M_{-1}^(3/4)      (vì -1/2 = 1/4*1 + 3/4*(-1)),   hoặc   omega^2 ~ sqrt(M_0 M_{-1}).

Script kiểm các khẳng định trên cho Duffing, x + x^7, x + x^3 + x^5, con lắc và Duffing-điều hòa, rồi đo sai số.
Chạy:  python power_mean_view.py   (đầu ra: power_mean_view.txt)
"""
import sys

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import j1

GL_X, GL_W = np.polynomial.legendre.leggauss(40)


def chord_stiffness(g, A):
    """K(theta) = 2/(A^2 (1-s^2)) int_{As}^{A} g = (1/(A(1+s))) int_s^1 g(Av) dv * 2/(1-s) ... ổn định số."""
    def K(theta):
        s = np.sin(theta)
        v = s + (1 - s) * (GL_X + 1) / 2                     # Gauss-Legendre trên [s, 1]
        integral = (1 - s) / 2 * np.dot(GL_W, g(A * v))       # = int_s^1 g(Av) dv
        return 2 * integral / (A * (1 - s * s)) if 1 - s > 1e-13 else g(A) / A
    return K


def E(f):
    return (2 / np.pi) * quad(f, 0, np.pi / 2, limit=400, epsabs=1e-15, epsrel=1e-13)[0]


def power_means(g, A):
    K = chord_stiffness(g, A)
    return {
        "M1": E(K),
        "M0": np.exp(E(lambda t: np.log(K(t)))),
        "Mh": E(lambda t: K(t) ** -0.5) ** -2,
        "Mm1": 1 / E(lambda t: 1 / K(t)),
    }


def omega_ode(g, A):
    rhs = lambda t, y: [y[1], -g(y[0])]
    ev = lambda t, y: y[1]
    ev.direction, ev.terminal = 1.0, True
    sol = solve_ivp(rhs, [0, 1e4], [A, 0.0], method="DOP853", rtol=1e-12, atol=1e-14, events=ev)
    return np.pi / sol.t_events[0][0]


def classical_time_average(g, A):
    """(2/A) <g(A cos t) cos t> trên một chu kỳ: công thức tuyến tính hóa cổ điển."""
    return (2 / A) * (1 / np.pi) * quad(lambda t: g(A * np.cos(t)) * np.cos(t), 0, np.pi, limit=200)[0]


CASES = [
    ("Duffing x+x^3", lambda x: x + x ** 3, [0.1, 0.5, 1, 2, 5, 10, 100, 1000]),
    ("x+|x|^7 sgn x", lambda x: x + np.abs(x) ** 7 * np.sign(x), [0.5, 1, 1.5, 2, 3, 10]),
    ("x+x^3+x^5", lambda x: x + x ** 3 + x ** 5, [0.5, 1, 2, 5, 20]),
    ("con lắc sin x", np.sin, [0.5, 1, 2, 2.5, 2.8, 3.0, 3.1]),
    ("Duffing-điều hòa x^3/(1+x^2)", lambda x: x ** 3 / (1 + x * x), [0.1, 0.5, 1, 2, 5, 50]),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    out = []
    P = out.append
    pct = lambda m, w: 100 * (np.sqrt(m) / w - 1)

    P("Sai số tần số (%) so với tích phân ODE. Cột 'kiểm' là sai số của M_{-1/2} (phải ~0).")
    P("%-30s %7s | %8s | %8s %8s %8s | %9s %9s" % ("hệ", "A", "kiểm", "M1 (cổ)", "M0", "M-1", "M1^¼M-1^¾",
                                                   "√(M0M-1)"))
    worst = {}
    for name, g, As in CASES:
        for A in As:
            m = power_means(g, A)
            w = omega_ode(g, A)
            chk = pct(m["Mh"], w)
            assert abs(chk) < 1e-7, (name, A, chk)
            assert abs(m["M1"] - classical_time_average(g, A)) < 1e-10 * max(1, m["M1"]), (name, A)
            assert m["Mm1"] < w * w < m["M0"] < m["M1"], (name, A)
            d1 = m["M1"] ** 0.25 * m["Mm1"] ** 0.75
            d2 = np.sqrt(m["M0"] * m["Mm1"])
            row = [pct(m["M1"], w), pct(m["M0"], w), pct(m["Mm1"], w), pct(d1, w), pct(d2, w)]
            P("%-30s %7g | %+8.1e | %+8.3f %+8.3f %+8.3f | %+9.4f %+9.4f" % (name, A, chk, *row))
            worst[name] = max(worst.get(name, 0), abs(row[3]))
    P("")
    P("Đã kiểm: M_{-1/2} = tần số đúng; M_1 = công thức cổ điển (2/A)<g(A cos) cos>; M_{-1} < omega^2 < M_0 < M_1.")

    # Duffing: mọi trung bình có dạng đóng (bước đầu AGM), nên công thức kép là dạng đóng
    P("")
    P("Duffing x'' + x + q x^3 = 0 (A = 1): M_1 = 1 + 3q/4, M_{-1} = sqrt((1+q/2)(1+q)),")
    P("  công thức kép  omega^2 ~ (1+3q/4)^(1/4) [(1+q/2)(1+q)]^(3/8).")
    q = np.logspace(-4, 8, 3000)
    from scipy.special import ellipk
    we = np.pi * np.sqrt(1 + q) / (2 * ellipk(q / (2 * (1 + q))))
    wd = ((1 + 0.75 * q) ** 0.25 * ((1 + q / 2) * (1 + q)) ** 0.375) ** 0.5
    e = 100 * (wd / we - 1)
    P("  sai số trên 1e-4 <= q <= 1e8: %+.4f%% .. %+.4f%%; giới hạn q -> inf: %+.4f%%"
      % (e.min(), e.max(), 100 * ((0.75 ** 0.25 * 0.5 ** 0.375) ** 0.5 / np.sqrt(0.7177700110) - 1)))
    # con lắc: M_1 = 2 J1(A)/A
    P("Con lắc: M_1 = 2 J_1(A)/A (kiểm: %.2e)." % abs(power_means(np.sin, 2.0)["M1"] - 2 * j1(2.0) / 2.0))

    text = "\n".join(out)
    print(text)
    open("power_mean_view.txt", "w", encoding="utf-8").write(text + "\n")


if __name__ == "__main__":
    main()
