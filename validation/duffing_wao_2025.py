"""Toán tử trung bình trọng số cục bộ (WLAO / WAO) cho dao động tự do Duffing.

Phương trình:  x'' + w0^2 x + gamma x^(2n+1) = 0,  x(0) = a, x'(0) = 0,  nghiệm thử x = a cos(2 pi tau).

Toán tử (Nguyen A.T., Anh N.D., Vietnam J. Mech. 42(3) (2020) 341-354, Eq. (6);
         Nguyen A.T., Linh N.N., Thang N.C., Ngoc N.A., Vinh L.Q., Anh N.D.,
         Int. J. Non-Linear Mech. 173 (2025) 105057, Eqs. (7)-(12), (23)-(24)):

    <f, h> = <f> + p B(f, h),
    B(f, h) = int_0^h f + int_{1-h}^1 f - 2h int_0^1 f,

  * bài 2020 (WLAO): p cố định (0.25, 0.5, 1), w = min_h w(p, h)                  -- Eq. (33)-(34)
  * bài 2025 (WEL):  p(h) = +-sqrt(1 - h^2),  w = min_{h, dấu} w(h)                  -- Eq. (23), (38)-(39)

Tần số (Eq. (35) bài 2025, Eq. (33) bài 2020 -- hai cách dẫn cho cùng một biểu thức):
    w^2(h) = <g(a cos 2 pi tau) cos 2 pi tau, h> / (a <cos^2 2 pi tau, h>).

Với Duffing (2n+1):  w^2(h) = w0^2 + gamma a^(2n) R_n(h),  R_n = <cos^(2n+2), h> / <cos^2, h>.
R_n không phụ thuộc a, gamma nên min theo h cho một hằng số lam_n duy nhất:
    w^2 = w0^2 + lam_n gamma a^(2n).

Script kiểm lại các bảng số trong hai bài báo rồi so sánh với các phương pháp ở
duffing_weighted_averaging.py (cổ điển 3/4, trọng số Anh 2017 s = 2, nghiệm chính xác).

Chạy:  python duffing_wao_2025.py
"""
import sys
from math import comb

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import ellipk

from duffing_weighted_averaging import lam_weighted, omega_exact_elliptic

H_GRID = np.linspace(0.0, 1.0, 2001)


# ---------- Toán tử ----------
def int_cos_pow(m, h):
    """int_0^h cos^m(2 pi tau) d tau, giải tích qua khai triển Fourier (m chẵn)."""
    total = comb(m, m // 2) / 2**m * h
    for j in range(m // 2):
        k = m - 2 * j
        total += 2 * comb(m, j) / 2**m * np.sin(2 * np.pi * k * h) / (2 * np.pi * k)
    return total


def B_cos_pow(m, h):
    """B(cos^m, h); với cos^m chẵn và tuần hoàn chu kỳ 1/2, int_{1-h}^1 = int_0^h."""
    return 2 * int_cos_pow(m, h) - 2 * h * int_cos_pow(m, 1.0)


def wao_cos_pow(m, h, p):
    """<cos^m(2 pi tau), h> = <cos^m> + p B(cos^m, h)."""
    return int_cos_pow(m, 1.0) + p * B_cos_pow(m, h)


def wao_numeric(f, h, p):
    """Kiểm chéo: định nghĩa gốc Eq. (7) với q = 1 - 2ph, tích phân số."""
    g = lambda t: f(t)
    full = quad(g, 0, 1, limit=200)[0]
    loc = quad(g, 0, h, limit=200)[0] + quad(g, 1 - h, 1, limit=200)[0]
    return (1 - 2 * p * h) * full + p * loc


def R_ratio(n, h, p):
    return wao_cos_pow(2 * n + 2, h, p) / wao_cos_pow(2, h, p)


def lam_min(n, p_of_h):
    """min_h R_n(h, p(h)), lưới mịn rồi tinh chỉnh bằng Brent. Trả về (lam, h*)."""
    vals = np.array([R_ratio(n, h, p_of_h(h)) for h in H_GRID])
    i = int(np.argmin(vals))
    lo, hi = H_GRID[max(i - 1, 0)], H_GRID[min(i + 1, len(H_GRID) - 1)]
    if hi - lo < 1e-12:
        return vals[i], H_GRID[i]
    res = minimize_scalar(lambda h: R_ratio(n, h, p_of_h(h)), bounds=(lo, hi), method="bounded",
                          options={"xatol": 1e-12})
    return (res.fun, res.x) if res.fun < vals[i] else (vals[i], H_GRID[i])


P_PLUS = lambda h: np.sqrt(max(1 - h * h, 0.0))
P_MINUS = lambda h: -np.sqrt(max(1 - h * h, 0.0))


def lam_wel_2025(n):
    """Bài 2025: min theo h trên cả hai nhánh p = +-sqrt(1-h^2). Trả về (lam, h*, nhánh)."""
    lp, hp = lam_min(n, P_PLUS)
    lm, hm = lam_min(n, P_MINUS)
    return (lp, hp, "+") if lp < lm else (lm, hm, "-")


def lam_wlao_2020(n, p):
    lam, h = lam_min(n, lambda _h: p)
    return lam, h


# ---------- Nghiệm chính xác Duffing bậc 2n+1 ----------
def omega_exact_general(a, gamma, n, w0=1.0):
    """Eq. (50) bài 2025, đổi biến x = a sin(phi):
    T = 4 int_0^{pi/2} a cos(phi) dphi / sqrt(2 (U(a) - U(x)))."""
    def U(x):
        return w0**2 * x**2 / 2 + gamma * x ** (2 * n + 2) / (2 * n + 2)

    def integrand(phi):
        x = a * np.sin(phi)
        # (U(a)-U(x)) / (a^2 cos^2 phi), viết dạng không kỳ dị
        s2 = np.sin(phi) ** 2
        geo = sum(s2**k for k in range(n + 1))          # (1 - s^(2n+2)) / (1 - s^2)
        val = w0**2 / 2 + gamma * a ** (2 * n) * geo / (2 * n + 2)
        return 1.0 / np.sqrt(2 * val)
    T = 4 * quad(integrand, 0, np.pi / 2, epsabs=1e-14, epsrel=1e-13)[0]
    return 2 * np.pi / T


def omega_lin(a, gamma, n, lam, w0=1.0):
    return np.sqrt(w0**2 + lam * gamma * a ** (2 * n))


def main():
    sys.stdout.reconfigure(encoding="utf-8")

    # 1. Kiểm tra toán tử: giải tích vs tích phân số, và vs Bảng 1 bài 2025
    for m in (2, 4, 6):
        for h in (0.07, 0.3, 0.42):
            for p in (0.5, -0.8):
                an = wao_cos_pow(m, h, p)
                nu = wao_numeric(lambda t: np.cos(2 * np.pi * t) ** m, h, p)
                assert abs(an - nu) < 1e-10, (m, h, p, an, nu)
    for h in (0.1, 0.37):
        for sgn in (1, -1):
            p = sgn * np.sqrt(1 - h * h)
            t2 = (2 * np.pi + p * np.sin(4 * np.pi * h)) / (4 * np.pi)
            t4 = (12 * np.pi + p * (np.sin(8 * np.pi * h) + 8 * np.sin(4 * np.pi * h))) / (32 * np.pi)
            assert abs(wao_cos_pow(2, h, p) - t2) < 1e-13
            assert abs(wao_cos_pow(4, h, p) - t4) < 1e-13
    # Tính chất (iii)/(v): B = 0 tại h = 0, 0.5, 1
    for h in (0.0, 0.5, 1.0):
        assert abs(B_cos_pow(4, h)) < 1e-13
    print("Kiểm tra toán tử: giải tích = tích phân số = Bảng 1 (IJNLM 2025); B(.,0)=B(.,0.5)=B(.,1)=0.\n")

    # 2. Tái lập Bảng 7 VJMech 2020: x'' + x^(2n+1) = 0, A = 1 (w/A^n)
    print("Tái lập Bảng 7 VJMech 2020 (x'' + x^(2n+1) = 0):  w/A^n = sqrt(lam_n)")
    print("%3s %8s | %8s %8s %8s %8s | %s" % ("n", "exact", "p=0", "p=0.25", "p=0.5", "p=1", "h* (p=0.25,0.5,1)"))
    for n in range(1, 6):
        we = omega_exact_general(1.0, 1.0, n, w0=0.0)
        row = [np.sqrt(int_cos_pow(2 * n + 2, 1) / 0.5)]
        hs = []
        for p in (0.25, 0.5, 1.0):
            lam, h = lam_wlao_2020(n, p)
            row.append(np.sqrt(lam))
            hs.append(h)
        print("%3d %8.4f | %8.4f %8.4f %8.4f %8.4f | %.3f %.3f %.3f" % (n, we, *row, *hs))
    print("  Bài báo: n=1: 0.8472 | 0.8660 0.8583 0.8498 0.8306 ; n=5: 0.5772 | 0.6717 0.6557 0.6383 0.5983\n")

    # 3. Tái lập Bảng 4-5 IJNLM 2025
    print("Tái lập Bảng 4-5 IJNLM 2025 (w0 = 1):  lam_n = min_{h,+-} R_n")
    lams = {}
    for n in range(1, 6):
        lam, h, br = lam_wel_2025(n)
        lams[n] = lam
        print("  n=%d: lam = %.5f, h* = %.3f, nhánh p = %ssqrt(1-h^2)" % (n, lam, h, br))
    print("%3s %5s %5s | %10s | %10s %7s" % ("n", "gam", "a", "w_exact", "w_WEL", "err%"))
    for gamma, a in ((1, 1), (10, 1), (10, 10)):
        for n in range(1, 6):
            we = omega_exact_general(a, gamma, n)
            ww = omega_lin(a, gamma, n, lams[n])
            print("%3d %5g %5g | %10.4g | %10.4g %+7.2f" % (n, gamma, a, we, ww, 100 * (ww - we) / we))
    print("  Bài báo (n=1): 1.300 (-1.34%), 2.811 (-1.94%), 26.29 (-1.94%); "
          "(n=5, gam=a=10): 189317 (+3.73%)\n")

    # 4. So sánh cho Duffing bậc 3 (n = 1) 
    lam_c = 0.75
    lam_s2 = lam_weighted(2.0)
    lam_2025 = lams[1]
    lam_2020 = {p: lam_wlao_2020(1, p)[0] for p in (0.25, 0.5, 1.0)}
    lam_inf = (np.pi / (2 * ellipk(0.5))) ** 2
    methods = [("cổ điển 3/4", lam_c), ("Anh 2017 s=2", lam_s2)]
    methods += [("WLAO 2020 p=%g" % p, l) for p, l in lam_2020.items()]
    methods += [("WEL 2025", lam_2025)]
    print("Duffing bậc 3: w^2 = w0^2 + lam beta A^2   (giới hạn chính xác lam_inf = %.4f)" % lam_inf)
    for name, lam in methods:
        print("  %-16s lam = %.4f" % (name, lam))
    print()

    q = np.logspace(-3, 6, 4000)                      # beta A^2 / w0^2
    we = np.array([omega_exact_elliptic(np.sqrt(v), 1.0) for v in q])
    print("Sai số w trên 1e-3 <= beta A^2/w0^2 <= 1e6:")
    print("  %-16s %10s %10s %12s" % ("phương pháp", "min err%", "max err%", "max |err|%"))
    errs = {}
    for name, lam in methods:
        e = 100 * (np.sqrt(1 + lam * q) - we) / we
        errs[name] = e
        print("  %-16s %+10.3f %+10.3f %12.3f" % (name, e.min(), e.max(), np.abs(e).max()))
    print()

    # Phương pháp nào chính xác nhất theo từng vùng
    names = [m[0] for m in methods]
    E = np.abs(np.array([errs[nm] for nm in names]))
    best = np.argmin(E, axis=0)
    print("Phương pháp có |sai số| nhỏ nhất theo vùng beta A^2:")
    start = 0
    for i in range(1, len(q) + 1):
        if i == len(q) or best[i] != best[start]:
            print("  %9.3g .. %-9.3g : %s" % (q[start], q[i - 1], names[best[start]]))
            start = i
    print()

    print("Bảng chi tiết (beta = 1):")
    hdr = "%6s | %9s | " % ("A", "w_exact") + " | ".join("%16s" % nm for nm in names)
    print(hdr)
    print("-" * len(hdr))
    for A in (0.1, 0.5, 1.0, 2.0, 5.0, 10.0):
        w_e = omega_exact_elliptic(A, 1.0)
        cells = []
        for _, lam in methods:
            w = np.sqrt(1 + lam * A * A)
            cells.append("%8.5f %+6.2f%%" % (w, 100 * (w - w_e) / w_e))
        print("%6.1f | %9.5f | " % (A, w_e) + " | ".join(cells))


if __name__ == "__main__":
    main()
