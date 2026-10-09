"""Tần số dao động tự do của Duffing  x'' + w0^2 x + beta x^3 = 0,  x(0)=A, x'(0)=0.

So sánh:
  1. Nghiệm chính xác: tích phân elliptic loại 1, kiểm chéo bằng RK (DOP853, tol 1e-12).
  2. Tuyến tính hóa / cân bằng điều hòa cổ điển (Krylov-Bogoliubov): w^2 = w0^2 + (3/4) beta A^2.
  3. Tuyến tính hóa với trung bình trọng số của N.D. Anh (Anh 2015; Anh, Hai, Hieu,
     Lat. Am. J. Solids Struct. 14 (2017) 1723-1740):
        <f>_s = int_0^inf s^2 theta e^{-s theta} f(theta) dtheta,   s = 2 trong bài 2017.
  4. Tiêu chuẩn đối ngẫu (Anh-Hiếu-Linh 2012) áp với trung bình điều hòa cổ điển
     (cột tham khảo, do chúng tôi tự áp dụng cho bài toán tất định).

Chạy:  python duffing_weighted_averaging.py
"""
import sys
import numpy as np
from scipy.special import ellipk
from scipy.integrate import solve_ivp, quad

W0 = 1.0
BETAS = [0.1, 1.0, 10.0]
AMPS = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]


# ---------- Toán tử trung bình trọng số ----------
def wavg_cos(n, s):
    """<cos(n theta)>_s = s^2 (s^2 - n^2) / (s^2 + n^2)^2   (biến đổi Laplace)."""
    return s**2 * (s**2 - n**2) / (s**2 + n**2) ** 2


def wavg_cos_power(m, s):
    """<cos^m theta>_s qua khai triển Fourier của cos^m (m chẵn)."""
    from math import comb
    total = comb(m, m // 2) / 2**m                     # số hạng hằng
    for j in range(m // 2):                            # các số hạng cos((m-2j) theta)
        total += 2 * comb(m, j) / 2**m * wavg_cos(m - 2 * j, s)
    return total


def wavg_numeric(f, s):
    """Kiểm chéo bằng tích phân số trực tiếp."""
    return quad(lambda t: s**2 * t * np.exp(-s * t) * f(t), 0, np.inf, limit=400)[0]


def lam_weighted(s=2.0):
    """k_eq = w0^2 + lam * beta * A^2 với lam = <cos^4>_s / <cos^2>_s."""
    return wavg_cos_power(4, s) / wavg_cos_power(2, s)


def lam_dual_classical():
    """Tiêu chuẩn đối ngẫu: alpha = E[AB]/E[B^2] / (2 - mu), mu = E[AB]^2/(E[A^2]E[B^2]),
    A = beta x^3, B = x, trung bình điều hòa: <cos^2>=1/2, <cos^4>=3/8, <cos^6>=5/16."""
    c2, c4, c6 = 1 / 2, 3 / 8, 5 / 16
    mu = c4**2 / (c2 * c6)
    return (c4 / c2) / (2 - mu)


# ---------- Tần số ----------
def omega_lin(A, beta, lam):
    return np.sqrt(W0**2 + lam * beta * A**2)


def omega_exact_elliptic(A, beta):
    """w = pi sqrt(w0^2 + beta A^2) / (2 K(m)),  m = beta A^2 / (2 (w0^2 + beta A^2))."""
    a2 = W0**2 + beta * A**2
    m = beta * A**2 / (2 * a2)
    return np.pi * np.sqrt(a2) / (2 * ellipk(m))


def omega_exact_quadrature(A, beta):
    """T = 4 int_0^A dx / sqrt(2(U(A)-U(x))), đổi biến x = A sin(phi) để khử kỳ dị."""
    def integrand(phi):
        # U(A)-U(x) = (A^2 - x^2)[w0^2/2 + beta (A^2 + x^2)/4]
        x = A * np.sin(phi)
        return 1.0 / np.sqrt(2 * (W0**2 / 2 + beta * (A**2 + x**2) / 4))
    T = 4 * quad(integrand, 0, np.pi / 2, epsabs=1e-14, epsrel=1e-13)[0]
    return 2 * np.pi / T


def omega_exact_ode(A, beta):
    """Tích phân số đến lần đầu x' đổi dấu từ âm sang dương (nửa chu kỳ: x = -A)."""
    rhs = lambda t, y: [y[1], -W0**2 * y[0] - beta * y[0] ** 3]
    ev = lambda t, y: y[1]
    ev.direction = 1.0
    ev.terminal = True
    T_guess = 2 * np.pi / omega_lin(A, beta, 0.75)
    sol = solve_ivp(rhs, [0, 2 * T_guess], [A, 0.0], method="DOP853", rtol=1e-12, atol=1e-14,
                    events=ev)
    half_T = sol.t_events[0][0]
    return np.pi / half_T


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    s = 2.0
    lam_w = lam_weighted(s)
    lam_d = lam_dual_classical()

    # Kiểm tra toán tử trọng số bằng tích phân số
    for m in (2, 4, 6):
        an, nu = wavg_cos_power(m, s), wavg_numeric(lambda t: np.cos(t) ** m, s)
        assert abs(an - nu) < 1e-10, (m, an, nu)
    print("Toán tử trọng số s = 2:  <cos^2> = %.4f, <cos^4> = %.4f, <cos^6> = %.4f"
          % tuple(wavg_cos_power(m, s) for m in (2, 4, 6)))
    print("Hệ số k_eq = w0^2 + lam*beta*A^2:  cổ điển 0.7500 | trọng số (s=2) %.4f | "
          "đối ngẫu %.4f" % (lam_w, lam_d))
    lam_inf = (np.pi / (2 * ellipk(0.5))) ** 2
    print("Giới hạn chính xác khi beta A^2 -> inf: w^2/(beta A^2) -> %.4f\n" % lam_inf)

    hdr = ("%6s %6s | %10s | %10s %8s | %10s %8s | %10s %8s | %9s"
           % ("beta", "A", "w_exact", "w_classic", "err%", "w_weight", "err%", "w_dual", "err%",
              "ODE chk"))
    for beta in BETAS:
        print(hdr)
        print("-" * len(hdr))
        for A in AMPS:
            we = omega_exact_elliptic(A, beta)
            wq = omega_exact_quadrature(A, beta)
            wo = omega_exact_ode(A, beta)
            chk = max(abs(wq - we), abs(wo - we)) / we
            row = [beta, A, we]
            for lam in (0.75, lam_w, lam_d):
                w = omega_lin(A, beta, lam)
                row += [w, 100 * (w - we) / we]
            print("%6.1f %6.1f | %10.6f | %10.6f %+8.3f | %10.6f %+8.3f | %10.6f %+8.3f | %9.1e"
                  % (*row, chk))
        print()

    # Sai số lớn nhất trên toàn bộ dải (chỉ phụ thuộc beta A^2 / w0^2)
    q = np.logspace(-3, 6, 4000)
    we = np.array([omega_exact_elliptic(np.sqrt(v), 1.0) for v in q])
    for name, lam in (("cổ điển 3/4", 0.75), ("trọng số s=2", lam_w), ("đối ngẫu", lam_d)):
        err = 100 * (omega_lin(np.sqrt(q), 1.0, lam) - we) / we
        print("Sai số w trên 1e-3 <= beta A^2/w0^2 <= 1e6, %-13s: từ %+.3f%% đến %+.3f%%"
              % (name, err.min(), err.max()))


if __name__ == "__main__":
    main()
