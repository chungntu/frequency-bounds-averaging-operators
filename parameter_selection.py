"""Chọn tham số trọng số của các toán tử trung bình từ hai giới hạn đúng.

Với  x'' + w0^2 x + gamma |x|^alpha sgn x = 0, mọi toán tử trung bình cố định cho w^2 = w0^2 + lam gamma A^(alpha-1),
và hằng số đúng lam(q) giảm ngặt từ lam0 xuống lamInf (xem NOTES_vi.md). Script này:

  1. tìm hằng số minimax lam_mm(alpha): hằng số làm nhỏ nhất sai số tần số lớn nhất trên 1e-3 <= q <= 1e7;
  2. đo vị trí của lam_mm trong khoảng (lamInf, lam0) và thử quy tắc dạng đóng lam = lamInf + c (lam0 - lamInf);
  3. tìm tham số làm toán tử cho đúng lam_mm:
       s*  của Anh-Hai-Hieu 2017 (nhánh s > argmin, chứa s = 2),
       p*  của Nguyen-Anh 2020 (WLAO, p cố định),
       n*  của họ siêu ellipse |p|^n + h^n = 1 trong Nguyen et al. 2025 (bài báo dùng n = 2);
  4. đo độ nhạy của WEL 2025 theo n.

Chạy:  python parameter_selection.py   (khoảng 1-2 phút; đầu ra: parameter_selection.txt)
"""
import sys

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import gamma as G

from frequency_bounds import (ALPHAS, LocalAverager, _min_over_h, lam0, lam_anh2017, lam_exact, lam_inf,
                              lam_wlao2020)

Q = np.logspace(-3, 7, 400)
C_RULE = 0.16


def lam0_closed(alpha):
    """lam0 = 2 <|cos|^(alpha+1)> = 2 Gamma((alpha+2)/2) / (sqrt(pi) Gamma((alpha+3)/2))."""
    return 2 * G((alpha + 2) / 2) / (np.sqrt(np.pi) * G((alpha + 3) / 2))


def lam_inf_closed(alpha):
    """Bình phương tần số của x'' + |x|^alpha sgn x = 0, A = 1."""
    k = np.sqrt(np.pi * (alpha + 1) / 2) * G((3 + alpha) / (2 * (alpha + 1))) / G(1 / (alpha + 1))
    return k * k


def max_abs_err(lc, lt):
    return np.abs(100 * (np.sqrt((1 + lc * Q) / (1 + lt * Q)) - 1)).max()


def lam_wel_family(n, num, den):
    """Toán tử Nguyen et al. 2025 với p = +-(1 - h^n)^(1/n); bài báo dùng n = 2."""
    best = np.inf
    for sg in (1, -1):
        pf = lambda h: sg * max(1 - h ** n, 0.0) ** (1 / n)
        val, _ = _min_over_h(lambda h: num.avg(h, pf(h)) / den.avg(h, pf(h)))
        best = min(best, val)
    return best


def solve(f, lo, hi, **kw):
    try:
        return brentq(f, lo, hi, **kw)
    except ValueError:
        return np.nan


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    out = []
    P = out.append
    den = LocalAverager(2)

    for a in (0.3, 3.0, 8.0):
        assert abs(lam0_closed(a) - lam0(a)) < 1e-11 and abs(lam_inf_closed(a) - lam_inf(a)) < 1e-10
    P("Hai giới hạn có dạng đóng (hàm Gamma), khớp tích phân số tới 1e-10.\n")

    P("%5s %8s %8s %8s %6s %7s | %7s %7s | %6s %6s %6s | %s"
      % ("alpha", "lam0", "lamInf", "lam_mm", "vị trí", "err_mm", "rule", "err", "s*", "p*", "n*",
         "WEL lam(n=1, 2, 4, 8)"))
    rows = []
    for a in ALPHAS:
        lt = np.array([lam_exact(q, a) for q in Q])
        l0, li = lam0(a), lam_inf(a)
        r = minimize_scalar(lambda lc: max_abs_err(lc, lt), bounds=(li, l0), method="bounded",
                            options={"xatol": 1e-10})
        lmm, emm = r.x, r.fun
        frac = (lmm - li) / (l0 - li)
        lrule = li + C_RULE * (l0 - li)
        erule = max_abs_err(lrule, lt)
        s_star = solve(lambda s: lam_anh2017(a, s) - lmm, 1.4, 6.0)
        num = LocalAverager(a + 1)
        p_star = solve(lambda p: lam_wlao2020(a, p, num, den)[0] - lmm, 0.0, 2.0, xtol=1e-6)
        n_star = solve(lambda n: lam_wel_family(n, num, den) - lmm, 0.2, 40.0, xtol=1e-4)
        wel = [lam_wel_family(n, num, den) for n in (1, 2, 4, 8)]
        P("%5g %8.5f %8.5f %8.5f %6.3f %6.3f%% | %7.5f %6.3f%% | %6.3f %6.3f %6.2f | %s"
          % (a, l0, li, lmm, frac, emm, lrule, erule, s_star, p_star, n_star,
             " ".join("%.4f" % w for w in wel)))
        rows.append((a, frac, emm, erule, s_star, p_star, n_star, wel))

    fr = np.array([r[1] for r in rows])
    P("")
    P("Vị trí của lam_mm trong khoảng (lamInf, lam0): %.3f đến %.3f." % (fr.min(), fr.max()))
    P("Quy tắc lam = lamInf + %.2f (lam0 - lamInf): sai số lớn nhất vượt minimax nhiều nhất %.3f điểm %%."
      % (C_RULE, max(r[3] - r[2] for r in rows)))
    spread = max((max(r[7]) - min(r[7])) / min(r[7]) for r in rows)
    P("WEL 2025: đổi n từ 1 đến 8 làm lam thay đổi nhiều nhất %.2f%%; tham số n gần như không điều chỉnh được lam."
      % (100 * spread))
    r3 = [r for r in rows if r[0] == 3.0][0]
    P("Duffing (alpha = 3): s* = %.3f (bài 2017 dùng s = 2), p* = %.3f (bài 2020 thử 0,25; 0,5; 1)."
      % (r3[4], r3[5]))

    text = "\n".join(out)
    print(text)
    open("parameter_selection.txt", "w", encoding="utf-8").write(text + "\n")


if __name__ == "__main__":
    main()
