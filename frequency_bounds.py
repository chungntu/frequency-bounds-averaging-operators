"""Chặn hai phía cho tần số dao động tự do  x'' + w0^2 x + gamma |x|^alpha sgn x = 0,  x(0)=A, x'(0)=0.

Viết w^2 = w0^2 + lam(q) gamma A^(alpha-1),  q = gamma A^(alpha-1) / w0^2  (lấy w0 = 1).

Biểu diễn chính xác (đổi biến x = A sin(theta), xem GHI_CHU.md):
    1/w = E[(1 + q H)^(-1/2)],   theta ~ đều trên [0, pi/2],
    H(theta) = h(sin theta),  h(s) = 2 (1 - s^(alpha+1)) / ((alpha+1)(1 - s^2)).
Do đó w^2 = M_{-1/2}(1 + q H) (trung bình lũy thừa bậc -1/2), lõm theo q, và
    lam0   = E[H]                 (giới hạn biên độ nhỏ; trùng tuyến tính hóa cổ điển),
    lamInf = M_{-1/2}(H) = (E[H^(-1/2)])^(-2)   (giới hạn biên độ lớn; dao động thuần lũy thừa),
    lam'(0) = -(3/4) Var(H).

Script:
  1. kiểm biểu diễn trên với tích phân elliptic (alpha = 3) và tích phân ODE (alpha = 0.5, 2.5, 7);
  2. kiểm lam0 = 2 <|cos|^(alpha+1)>, lamInf = tần số dao động thuần lũy thừa, lam'(0);
  3. kiểm số tính đơn điệu của lam(q) trên lưới;
  4. xếp chỗ bốn toán tử trung bình (cổ điển; Anh-Hai-Hieu 2017 s=2; Nguyen-Anh 2020 WLAO p=0.25, 0.5, 1;
     Nguyen et al. 2025 WEL) theo alpha, kèm hằng số hữu tỷ qua hai giới hạn và độ dốc ban đầu;
  5. đối chiếu với Bảng 3 và Bảng 5 của bài IJNLM 2025.

Chạy:  python frequency_bounds.py   (đầu ra: output.txt, ranking.csv)
"""
import csv
import sys

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.special import ellipk, gamma as G

ALPHAS = [0.1, 0.25, 0.5, 0.75, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 9.0, 10.0, 11.0]
Q_GRID = np.logspace(-3, 6, 400)


# ---------- Biểu diễn chính xác ----------
def H(theta, alpha):
    s = np.sin(theta)
    s2 = s * s
    with np.errstate(invalid="ignore", divide="ignore"):
        v = 2 * (1 - s ** (alpha + 1)) / ((alpha + 1) * (1 - s2))
    return np.where(np.abs(1 - s2) < 1e-12, 1.0, v)       # giới hạn s -> 1 bằng 1


def mean_theta(f):
    return (2 / np.pi) * quad(f, 0, np.pi / 2, limit=400, epsabs=1e-14, epsrel=1e-12)[0]


def omega2_exact(q, alpha):
    m = mean_theta(lambda t: (1 + q * H(t, alpha)) ** -0.5)
    return m ** -2


def lam_exact(q, alpha):
    return (omega2_exact(q, alpha) - 1) / q


def lam0(alpha):
    return mean_theta(lambda t: H(t, alpha))


def lam_inf(alpha):
    return mean_theta(lambda t: H(t, alpha) ** -0.5) ** -2


def var_H(alpha):
    m1 = lam0(alpha)
    return mean_theta(lambda t: H(t, alpha) ** 2) - m1 * m1


# ---------- Kiểm chéo độc lập ----------
def omega_ode(q, alpha):
    """Tích phân trực tiếp x'' + x + q |x|^alpha sgn x = 0, x(0)=1: thời điểm x đạt -1."""
    rhs = lambda t, y: [y[1], -y[0] - q * np.abs(y[0]) ** alpha * np.sign(y[0])]
    ev = lambda t, y: y[1]
    ev.direction, ev.terminal = 1.0, True
    T0 = 2 * np.pi / np.sqrt(1 + lam0(alpha) * q)
    sol = solve_ivp(rhs, [0, 3 * T0], [1.0, 0.0], method="DOP853", rtol=1e-12, atol=1e-14, events=ev)
    return np.pi / sol.t_events[0][0]


def omega_ode_mixed(q3, q5):
    rhs = lambda t, y: [y[1], -y[0] - q3 * y[0] ** 3 - q5 * y[0] ** 5]
    ev = lambda t, y: y[1]
    ev.direction, ev.terminal = 1.0, True
    T0 = 2 * np.pi / np.sqrt(1 + 0.75 * q3 + 0.625 * q5)
    sol = solve_ivp(rhs, [0, 3 * T0], [1.0, 0.0], method="DOP853", rtol=1e-12, atol=1e-14, events=ev)
    return np.pi / sol.t_events[0][0]


def cos_power_mean(alpha):
    """<|cos|^(alpha+1)> trung bình đều trên một chu kỳ."""
    return G((alpha + 2) / 2) / (np.sqrt(np.pi) * G((alpha + 3) / 2))


def lam_inf_pure(alpha):
    """Tần số chính xác dao động thuần x'' + |x|^alpha sgn x = 0, A = 1, alpha >= 1 (Eq. (43) bài IJNLM 2025)."""
    k =np.sqrt(np.pi * (alpha + 1) / 2) * G((3 + alpha) / (2 * (alpha + 1))) \
        / G(1 / (alpha + 1))
    return k * k


# ---------- Bốn toán tử ----------
def lam_anh2017(alpha, s=2.0):
    """<f>_s = int_0^inf s^2 t e^{-s t} f(t) dt; dùng chu kỳ pi của |cos|^(alpha+1) và cos^2."""
    def avg(f):
        # tổng theo chu kỳ: int_0^inf w(t) f(t) dt = sum_k int_0^pi w(t + k pi) f(t) dt
        r = np.exp(-s * np.pi)
        # sum_k (t + k pi) e^{-s(t+k pi)} = e^{-s t} [t/(1-r) + pi r/(1-r)^2]
        wsum = lambda t: s * s * np.exp(-s * t) * (t / (1 - r) + np.pi * r / (1 - r) ** 2)
        return quad(lambda t: wsum(t) * f(t), 0, np.pi, limit=400, epsabs=1e-14, epsrel=1e-12)[0]
    num = avg(lambda t: np.abs(np.cos(t)) ** (alpha + 1))
    den = avg(lambda t: np.cos(t) ** 2)
    return num / den


class LocalAverager:
    """<f,h> = <f> + p B(f,h) cho f = |cos 2 pi tau|^m, tính trên lưới tau rất mịn."""
    N = 400001

    def __init__(self, m):
        self.tau = np.linspace(0, 1, self.N)
        f = np.abs(np.cos(2 * np.pi * self.tau)) ** m
        self.F = cumulative_trapezoid(f, self.tau, initial=0.0)
        self.full = self.F[-1]

    def I(self, h):
        return np.interp(h, self.tau, self.F)

    def avg(self, h, p):
        # |cos 2 pi tau| đối xứng quanh tau = 1/2 nên int_{1-h}^1 = int_0^h
        return self.full + p * (2 * self.I(h) - 2 * h * self.full)


def _min_over_h(R):
    hs = np.linspace(0, 1, 4001)
    v = np.array([R(h) for h in hs])
    i = int(np.argmin(v))
    lo, hi = hs[max(i - 1, 0)], hs[min(i + 1, len(hs) - 1)]
    res = minimize_scalar(R, bounds=(lo, hi), method="bounded", options={"xatol": 1e-10})
    return (res.fun, res.x) if res.fun < v[i] else (v[i], hs[i])


def lam_wlao2020(alpha, p, num, den):
    return _min_over_h(lambda h: num.avg(h, p) / den.avg(h, p))


def lam_wel2025(alpha, num, den):
    best = None
    for sgn in (1, -1):
        val, h = _min_over_h(lambda h: num.avg(h, sgn * np.sqrt(max(1 - h * h, 0))) /
                             den.avg(h, sgn * np.sqrt(max(1 - h * h, 0))))
        if best is None or val < best[0]:
            best = (val, h, "+" if sgn > 0 else "-")
    return best


# ---------- Phân loại ----------
def classify(lam_c, l0, li, alpha):
    if lam_c >= l0 - 1e-12:
        return "trên/bằng lam0", "luôn dương", np.nan
    if lam_c <= li + 1e-12:
        return "dưới lamInf", "luôn âm", np.nan
    qs = brentq(lambda lq: lam_exact(10 ** lq, alpha) - lam_c, -8, 12, xtol=1e-10)
    return "trong khoảng", "đổi dấu một lần", 10 ** qs


def max_err(lam_fun, alpha, w2_exact):
    w = np.sqrt(1 + np.array([lam_fun(q) for q in Q_GRID]) * Q_GRID)
    e = 100 * (w - np.sqrt(w2_exact)) / np.sqrt(w2_exact)
    return e.min(), e.max()


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    out = []
    P = lambda *a: out.append(" ".join(str(x) for x in a))

    # 1. Kiểm biểu diễn chính xác
    for A, beta in ((0.5, 1.0), (2.0, 1.0), (10.0, 10.0)):
        q = beta * A * A
        a2 = 1 + q
        m = q / (2 * a2)
        we = np.pi * np.sqrt(a2) / (2 * ellipk(m))
        assert abs(np.sqrt(omega2_exact(q, 3.0)) - we) / we < 1e-11
    worst = 0.0
    for alpha in (0.5, 2.5, 7.0):
        for q in (0.1, 3.0, 50.0):
            wr = np.sqrt(omega2_exact(q, alpha))
            worst = max(worst, abs(wr - omega_ode(q, alpha)) / wr)
    assert worst < 1e-9, worst
    P("1. Biểu diễn 1/w = E[(1+qH)^(-1/2)]: khớp elliptic (alpha=3) tới 1e-11, khớp ODE "
      "(alpha=0.5, 2.5, 7) tới %.1e." % worst)

    # 2. Kiểm hai giới hạn và độ dốc
    for alpha in (0.25, 1.5, 3.0, 7.0, 10.0):
        assert abs(lam0(alpha) - 2 * cos_power_mean(alpha)) < 1e-11, alpha
        if alpha >= 1:
            assert abs(lam_inf(alpha) - lam_inf_pure(alpha)) < 1e-10, alpha
        qq = 1e-5
        slope = (lam_exact(qq, alpha) - lam0(alpha)) / qq
        assert abs(slope + 0.75 * var_H(alpha)) < 2e-4 * max(1, var_H(alpha)), (alpha, slope)
    assert abs(var_H(3.0) - 1 / 32) < 1e-12
    P("2. lam0 = 2<|cos|^(alpha+1)> (tuyến tính hóa cổ điển); lamInf = tần số dao động thuần lũy thừa "
      "(công thức Gamma); lam'(0) = -(3/4)Var(H). Duffing: lam0 = 3/4, lamInf = %.6f, "
      "lam'(0) = -3/128." % lam_inf(3.0))

    # 3. Đơn điệu trên lưới
    worst_inc = -np.inf
    for alpha in ALPHAS:
        lv = np.array([lam_exact(q, alpha) for q in Q_GRID])
        worst_inc = max(worst_inc, np.max(np.diff(lv)))
        assert np.all(np.diff(lv) < 0), alpha
        assert lam_inf(alpha) < lv.min() and lv.max() < lam0(alpha)
    P("3. lam(q) giảm ngặt trên 400 điểm 1e-3 <= q <= 1e6 với mọi alpha trong %s "
      "(bước tăng lớn nhất %.2e < 0) và nằm trong (lamInf, lam0)." % (ALPHAS, worst_inc))

    # 3b. Đa thức hỗn hợp x + q3 x^3 + q5 x^5 (Định lý 3) và lò xo mềm (Hệ quả 2)
    l03, li3, l05, li5 = lam0(3.0), lam_inf(3.0), lam0(5.0), lam_inf(5.0)
    n_ok = 0
    gap_lo, gap_hi = np.inf, np.inf
    for q3 in (0.0, 0.1, 1.0, 10.0, 100.0):
        for q5 in (1e-3, 0.1, 1.0, 10.0, 1e3):
            w2 = (mean_theta(lambda t: (1 + q3 * H(t, 3.0) + q5 * H(t, 5.0)) ** -0.5)) ** -2
            w2_ode = omega_ode_mixed(q3, q5) ** 2
            assert abs(w2 - w2_ode) / w2 < 1e-9
            lo, hi = 1 + li3 * q3 + li5 * q5, 1 + l03 * q3 + l05 * q5
            assert lo < w2 < hi, (q3, q5)
            gap_lo, gap_hi = min(gap_lo, w2 - lo), min(gap_hi, hi - w2)
            n_ok += 1
    P("3b. Đa thức hỗn hợp x + q3 x^3 + q5 x^5: chặn hai phía đúng ở %d/%d điểm (biểu diễn khớp ODE "
      "tới 1e-9)." % (n_ok, n_ok))
    soft = []
    for alpha, qs in ((3.0, (-0.2, -0.5, -0.9)), (0.5, (-0.3, -0.6))):
        for q in qs:
            we = np.sqrt(omega2_exact(q, alpha))
            wo = omega_ode(q, alpha)
            wc = np.sqrt(1 + lam0(alpha) * q)
            assert abs(we - wo) / we < 1e-9 and wc > we
            soft.append("alpha=%g q=%g: %+.3f%%" % (alpha, q, 100 * (wc - we) / we))
    P("3c. Lò xo mềm (q < 0, trước khi mất dao động): cổ điển vẫn đánh giá cao tần số: "
      + "; ".join(soft))
    P("")

    # 4. Xếp chỗ các toán tử
    den = LocalAverager(2)
    rows = []
    P("4. Xếp chỗ các toán tử. Sai số tần số trên 1e-3 <= q <= 1e6 (min/max, %).")
    P("   Quy tắc: lam_op >= lam0 -> luôn đánh giá cao; lam_op <= lamInf -> luôn đánh giá thấp;"
      " ở giữa -> đúng tại một q* duy nhất.")
    for alpha in ALPHAS:
        l0, li = lam0(alpha), lam_inf(alpha)
        w2 = np.array([omega2_exact(q, alpha) for q in Q_GRID])
        num = LocalAverager(alpha + 1)
        ops = [("Cổ điển", l0)]
        ops.append(("Anh 2017 s=2", lam_anh2017(alpha)))
        for p in (0.25, 0.5, 1.0):
            ops.append(("WLAO 2020 p=%g" % p, lam_wlao2020(alpha, p, num, den)[0]))
        lw, hw, br = lam_wel2025(alpha, num, den)
        ops.append(("WEL 2025", lw))
        b = 0.75 * var_H(alpha) / (l0 - li)
        P("")
        P("alpha = %-5g  lam0 = %.5f  lamInf = %.5f  (độ rộng khoảng %.2f%% của lam0)"
          % (alpha, l0, li, 100 * (l0 - li) / l0))
        P("   %-17s %9s  %-15s %-16s %10s %16s" % ("toán tử", "lam", "vị trí", "dấu sai số w", "q*",
                                                   "sai số w (%)"))
        for name, lc in ops:
            pos, sign, qs = classify(lc, l0, li, alpha)
            emin, emax = max_err(lambda q, lc=lc: lc, alpha, w2)
            P("   %-17s %9.5f  %-15s %-16s %10s %+7.3f..%+7.3f"
              % (name, lc, pos, sign, "-" if np.isnan(qs) else "%.3g" % qs, emin, emax))
            rows.append([alpha, name, lc, l0, li, pos, sign, qs, emin, emax])
        for name, lf in (("lamInf (cận)", lambda q: li),
                         ("Hữu tỷ 2 đầu", lambda q: (l0 + li * b * q) / (1 + b * q))):
            emin, emax = max_err(lf, alpha, w2)
            P("   %-17s %9s  %-15s %-16s %10s %+7.3f..%+7.3f"
              % (name, "%.5f" % li if name.startswith("lamInf") else "b=%.3f" % b, "", "", "",
                 emin, emax))
            rows.append([alpha, name, li if name.startswith("lamInf") else np.nan, l0, li, "", "",
                         np.nan, emin, emax])

    # 5. Đối chiếu bài IJNLM 2025
    P("")
    P("5. Đối chiếu IJNLM 2025.")
    P("   Bảng 3 (dao động thuần, sai số = sqrt(lam_op/lamInf) - 1):")
    for alpha, paper in ((0.5, (0.31, -0.81)), (3.0, (2.22, -1.93)), (10.0, (14.76, 2.87))):
        num = LocalAverager(alpha + 1)
        li = lam_inf(alpha)
        e_el = 100 * (np.sqrt(lam0(alpha) / li) - 1)
        e_wel = 100 * (np.sqrt(lam_wel2025(alpha, num, den)[0] / li) - 1)
        P("     alpha=%-4g EL %+6.2f%% (bài %+.2f)   WEL %+6.2f%% (bài %+.2f)"
          % (alpha, e_el, paper[0], e_wel, paper[1]))
    P("   Bảng 5, gamma = a = 10 (q = 10^(2n+1)), Anh 2017 s=2 (bài: 0.15 1.53 3.02 4.42 5.70 %):")
    for n in range(1, 6):
        alpha = 2 * n + 1
        q = 10.0 ** (2 * n + 1)
        we = np.sqrt(omega2_exact(q, alpha))
        wa = np.sqrt(1 + lam_anh2017(alpha) * q)
        P("     n=%d  w_exact=%.6g  w_AW=%.6g  sai số %+.2f%%" % (n, we, wa, 100 * (wa - we) / we))

    text = "\n".join(out)
    print(text)
    open("output.txt", "w", encoding="utf-8").write(text + "\n")
    with open("ranking.csv", "w", newline="", encoding="utf-8") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["alpha", "operator", "lambda", "lambda0", "lambda_inf", "position", "error_sign",
                      "q_star", "err_min_pct", "err_max_pct"])
        wtr.writerows(rows)


if __name__ == "__main__":
    main()
