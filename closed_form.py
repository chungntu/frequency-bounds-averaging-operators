"""Công thức tần số dạng hữu tỷ cho  x'' + w0^2 x + gamma |x|^alpha sgn x = 0, đúng ở cả hai đầu.

Với q = gamma A^(alpha-1) / w0^2 và H(theta) như trong frequency_bounds.py, hằng số đúng lam(q) có:
    q -> 0 :  lam = L0 + L1 q + L2 q^2 + ...,   L0 = mu1,  L1 = -(3/4)(mu2 - mu1^2),
              L2 = (5/8) mu3 - (9/8) mu1 mu2 + (1/2) mu1^3,           mu_k = E[H^k]
    q -> oo:  lam = Linf + Cinf / q + ...,       Linf = (E H^(-1/2))^(-2),  Cinf = Linf^(3/2) E[H^(-3/2)] - 1
(khai triển của omega^2 = (E[(1+qH)^(-1/2)])^(-2) quanh q = 0 và quanh q = oo).

Padé hai điểm [2/2]:
    lam(q) ~ (L0 + a1 q + Linf b2 q^2) / (1 + b1 q + b2 q^2)
khớp đúng 5 thông tin trên; a1, b1, b2 giải được dạng đóng từ chúng. Không có tham số khớp số nào.
    omega^2 = w0^2 + gamma A^(alpha-1) lam(q).
L0 và Linf có dạng Gamma; mu2, mu3 và E[H^(-3/2)] là tích phân một lớp của hàm sơ cấp, lập bảng một lần cho mỗi alpha.

Chạy:  python closed_form.py   (đầu ra: closed_form.txt, closed_form_coefficients.csv)
"""
import csv
import sys

import numpy as np
from scipy.special import ellipk

from frequency_bounds import H, lam_exact, lam_inf, mean_theta

ALPHA_TABLE = [0.1, 0.25, 0.5, 0.75, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 13.0, 15.0, 20.0]
Q = np.logspace(-5, 9, 700)


def coefficients(alpha):
    mu1, mu2, mu3 = (mean_theta(lambda t, k=k: H(t, alpha) ** k) for k in (1, 2, 3))
    L0 = mu1
    L1 = -0.75 * (mu2 - mu1 ** 2)
    L2 = 5 * mu3 / 8 - 9 * mu1 * mu2 / 8 + mu1 ** 3 / 2
    Li = lam_inf(alpha)
    Ci = Li ** 1.5 * mean_theta(lambda t: H(t, alpha) ** -1.5) - 1
    d = L0 - Li
    b1 = (L2 + d * L1 / Ci) / (-d * d / Ci - L1)
    a1 = L1 + L0 * b1
    b2 = (L1 + d * b1) / Ci
    return dict(L0=L0, L1=L1, L2=L2, Linf=Li, Cinf=Ci, a1=a1, b1=b1, b2=b2)


def lam_pade(q, c):
    return (c["L0"] + c["a1"] * q + c["Linf"] * c["b2"] * q * q) / (1 + c["b1"] * q + c["b2"] * q * q)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    out = []
    P = out.append

    # Kiểm các hệ số khai triển bằng sai phân của lam_exact
    for a in (0.5, 3.0, 7.0):
        c = coefficients(a)
        q = 1e-3
        series = c["L0"] + c["L1"] * q + c["L2"] * q * q
        assert abs(lam_exact(q, a) - series) < 5e-9, a
        qb = 1e6
        assert abs((lam_exact(qb, a) - c["Linf"]) * qb - c["Cinf"]) < 1e-4, a
    c3 = coefficients(3.0)
    assert abs(c3["L1"] + 3 / 128) < 1e-12 and abs(c3["L2"] - 9 / 512) < 1e-12
    P("Kiểm: khai triển tại q = 0 (tới q^2) và tại q = oo (tới 1/q) khớp lam_exact. Duffing: L1 = -3/128, L2 = 9/512.")
    P("")
    P("%6s %9s %10s %9s %9s %9s %9s %9s %9s | %s"
      % ("alpha", "L0", "L1", "L2", "Linf", "Cinf", "a1", "b1", "b2", "sai số tần số (%) trên 1e-5<=q<=1e9"))
    rows = []
    for a in ALPHA_TABLE:
        c = coefficients(a)
        assert c["b1"] > 0 and c["b2"] > 0
        lt = np.array([lam_exact(q, a) for q in Q])
        w = np.sqrt(1 + lt * Q)
        e = 100 * (np.sqrt(1 + lam_pade(Q, c) * Q) / w - 1)
        P("%6g %9.6f %+10.6f %+9.6f %9.6f %9.6f %+9.5f %9.5f %9.5f | %+.6f .. %+.6f"
          % (a, c["L0"], c["L1"], c["L2"], c["Linf"], c["Cinf"], c["a1"], c["b1"], c["b2"], e.min(), e.max()))
        rows.append([a] + [c[k] for k in ("L0", "L1", "L2", "Linf", "Cinf", "a1", "b1", "b2")] + [e.min(), e.max()])

    # So sánh với các công thức Duffing đã biết
    P("")
    P("Duffing (alpha = 3), sai số tần số lớn nhất trên mọi biên độ:")
    q = np.logspace(-5, 9, 4000)
    we = np.pi * np.sqrt(1 + q) / (2 * ellipk(q / (2 * (1 + q))))
    m = q / (2 * (1 + q))
    agm2 = np.sqrt(1 + q) * (1 + (1 - m) ** 0.25) ** 2 / 4
    cmp = [("cổ điển 3/4", np.sqrt(1 + 0.75 * q)),
           ("AGM hai bước (đã biết)", agm2),
           ("Padé hai điểm (công thức này)", np.sqrt(1 + lam_pade(q, c3) * q))]
    for name, w in cmp:
        P("  %-32s %.6f%%" % (name, np.abs(100 * (w / we - 1)).max()))
    P("  Beléndez et al. 2009 (cân bằng điều hòa hữu tỷ bậc 2, theo bài gốc): 0.0055%")

    text = "\n".join(out)
    print(text)
    open("closed_form.txt", "w", encoding="utf-8").write(text + "\n")
    with open("closed_form_coefficients.csv", "w", newline="", encoding="utf-8") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["alpha", "L0", "L1", "L2", "Linf", "Cinf", "a1", "b1", "b2", "err_min_pct", "err_max_pct"])
        wtr.writerows(rows)


if __name__ == "__main__":
    main()
