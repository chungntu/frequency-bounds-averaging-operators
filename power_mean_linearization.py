"""Tuyến tính hóa bằng trung bình lũy thừa: thay trung bình cộng bằng trung bình bậc r.

Khung tuyến tính hóa cổ điển với nghiệm thử x = A cos(t):
    k = <g(x) x> / <x^2>  =  trung bình cộng của độ cứng cát tuyến k_s = g(x)/x với trọng số x^2 = A^2 cos^2 t.
Với g = w0^2 x + gamma |x|^alpha sgn x:  k_s / w0^2 = 1 + q |cos t|^(alpha-1),  q = gamma A^(alpha-1)/w0^2.

Đề xuất: giữ nguyên nghiệm thử và trọng số, chỉ thay trung bình cộng bằng trung bình lũy thừa bậc r:
    omega^2 / w0^2 = M_r^{(cos^2)}(1 + q |cos t|^(alpha-1)) = ( E_w[(1 + q|cos|^(alpha-1))^r] )^(1/r),
với E_w là trung bình theo trọng số cos^2 t. r = 1 là tuyến tính hóa cổ điển. Khác mọi toán tử tuyến tính,
trung bình bậc r != 1 cho lam phụ thuộc biên độ.

Script:
  1. bậc r đúng ở hai đầu: r0 (khớp độ dốc tại q -> 0, r0 = 1 - 1.5 Var_theta(H) / Var_w(|cos|^(alpha-1)))
     và rInf (khớp lamInf khi q -> oo);
  2. bậc r cố định tốt nhất cho từng alpha và sai số của một bậc chung r = 0.44;
  3. so với hằng số minimax (mục 8) và cách cổ điển.

Chạy:  python power_mean_linearization.py   (đầu ra: power_mean_linearization.txt)
"""
import sys

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from frequency_bounds import ALPHAS, lam0, lam_inf, omega2_exact, var_H

_X, _WQ = np.polynomial.legendre.leggauss(400)
T = (_X + 1) * np.pi / 4
W = _WQ * np.pi / 4 * np.cos(T) ** 2
W = W / W.sum()
Q = np.logspace(-3, 7, 160)
R_UNIVERSAL = 0.44


def power_mean_w(values, r):
    if abs(r) < 1e-12:
        return np.exp(np.sum(W * np.log(values)))
    return np.sum(W * values ** r) ** (1 / r)


def omega2_power_mean(q, alpha, r):
    return power_mean_w(1 + q * np.abs(np.cos(T)) ** (alpha - 1), r)


def max_err(alpha, r, w2):
    w = np.sqrt(np.array([omega2_power_mean(q, alpha, r) for q in Q]))
    return np.abs(100 * (w / np.sqrt(w2) - 1)).max()


def max_err_const(lc, w2):
    return np.abs(100 * (np.sqrt((1 + lc * Q) / w2) - 1)).max()


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    out = []
    P = out.append

    # r = 1 phải trùng tuyến tính hóa cổ điển (lam0)
    for a in (0.5, 3.0, 7.0):
        assert abs((omega2_power_mean(2.0, a, 1.0) - 1) / 2.0 - lam0(a)) < 1e-10
    P("Kiểm: r = 1 trùng tuyến tính hóa cổ điển (lam0).")
    P("")
    P("%6s %7s %7s | %7s %8s | %8s | %9s %9s" % ("alpha", "r0", "rInf", "r*", "err(r*)", "err(0.44)",
                                                  "minimax λ", "cổ điển"))
    rows = []
    for a in ALPHAS:
        w2 = np.array([omega2_exact(q, a) for q in Q])
        X = np.abs(np.cos(T)) ** (a - 1)
        varw = np.sum(W * X ** 2) - np.sum(W * X) ** 2
        r0 = 1 - 1.5 * var_H(a) / varw
        rinf = brentq(lambda r: power_mean_w(X, r) - lam_inf(a), -0.99, 3.0)
        res = minimize_scalar(lambda r: max_err(a, r, w2), bounds=(-0.5, 1.2), method="bounded",
                              options={"xatol": 1e-4})
        e_u = max_err(a, R_UNIVERSAL, w2)
        rmm = minimize_scalar(lambda lc: max_err_const(lc, w2), bounds=(lam_inf(a), lam0(a)), method="bounded",
                              options={"xatol": 1e-10})
        e_cl = max_err_const(lam0(a), w2)
        P("%6g %7.4f %7.4f | %7.3f %7.3f%% | %7.3f%% | %8.3f%% %8.3f%%" % (a, r0, rinf, res.x, res.fun, e_u, rmm.fun,
                                                                       e_cl))
        rows.append((a, res.x, res.fun, e_u, rmm.fun))
    rs = np.array([r[1] for r in rows])
    P("")
    P("Bậc tối ưu r* nằm trong [%.3f, %.3f] với mọi alpha; bậc chung r = %.2f kém tối ưu nhiều nhất %.3f điểm %%."
      % (rs.min(), rs.max(), R_UNIVERSAL, max(r[3] - r[2] for r in rows)))
    P("Bậc chung r = %.2f tốt hơn hằng số minimax (tốt nhất trong mọi toán tử tuyến tính) ở %d/%d giá trị alpha."
      % (R_UNIVERSAL, sum(r[3] < r[4] for r in rows), len(rows)))

    text = "\n".join(out)
    print(text)
    open("power_mean_linearization.txt", "w", encoding="utf-8").write(text + "\n")


if __name__ == "__main__":
    main()
