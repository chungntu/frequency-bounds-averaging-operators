"""Giải thích vị trí của hằng số minimax trong khoảng (lamInf, lam0) (mục 8: khoảng 0.14-0.17).

Viết lam(q) = lamInf + D phi(q), D = lam0 - lamInf, phi giảm từ 1 về 0. Với một hằng số lam_c = lamInf + delta D,
sai số tần số khi D nhỏ là  e(q) ~ (D/2) (delta - phi(q)) u(q) / lamInf,  u = lamInf q / (1 + lamInf q) in (0,1).
Sai số dương lớn nhất đạt khi q -> oo (u -> 1): delta. Sai số âm lớn nhất: max_u (phi - delta) u.
Minimax cân bằng hai giá trị:   delta = max_u (phi(u) - delta) u.

Mô hình hai giới hạn (mục 5): phi = 1/(1 + b q), b = (3/4) Var(H) / D, đổi sang biến u:
    phi(u) = (1-u) / (1 - u + beta u),   beta = b / lamInf = (3/4) Var(H) / (lamInf (lam0 - lamInf)).
Khi beta = 1, phi = 1 - u và  delta = (1-delta)^2/4,  tức  delta = 3 - 2 sqrt(2) = 0.17157...

Chạy:  python minimax_explained.py   (đầu ra: minimax_explained.txt)
"""
import sys

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from frequency_bounds import ALPHAS, lam0, lam_exact, lam_inf, var_H

U = np.linspace(0, 1, 200001)[1:-1]
Q = np.logspace(-3, 7, 400)


def delta_model(beta):
    phi = (1 - U) / (1 - U + beta * U)
    return brentq(lambda d: d - np.max((phi - d) * U), 1e-9, 0.9)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    out = []
    P = out.append
    d1 = delta_model(1.0)
    assert abs(d1 - (3 - 2 * np.sqrt(2))) < 1e-9
    P("Mô hình với beta = 1: delta = %.10f = 3 - 2 sqrt(2)." % d1)
    P("")
    P("%6s %8s | %10s %10s %10s" % ("alpha", "beta", "đo được", "mô hình", "lệch"))
    worst = 0.0
    for a in ALPHAS:
        l0, li = lam0(a), lam_inf(a)
        beta = 0.75 * var_H(a) / (li * (l0 - li))
        lt = np.array([lam_exact(q, a) for q in Q])
        f = lambda lc: np.abs(np.sqrt((1 + lc * Q) / (1 + lt * Q)) - 1).max()
        lmm = minimize_scalar(f, bounds=(li, l0), method="bounded", options={"xatol": 1e-12}).x
        frac = (lmm - li) / (l0 - li)
        dm = delta_model(beta)
        worst = max(worst, abs(frac - dm))
        P("%6g %8.4f | %10.4f %10.4f %+10.4f" % (a, beta, frac, dm, dm - frac))
    P("")
    P("Lệch lớn nhất giữa vị trí minimax đo được và mô hình: %.4f." % worst)
    P("beta gần 1 với 0.5 <= alpha <= 5 nên delta gần 3 - 2 sqrt(2) = 0.1716; beta tăng khi alpha lớn nên delta giảm.")

    text = "\n".join(out)
    print(text)
    open("minimax_explained.txt", "w", encoding="utf-8").write(text + "\n")


if __name__ == "__main__":
    main()
