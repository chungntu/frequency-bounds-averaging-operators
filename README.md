# Two-sided frequency bounds for nonlinear free oscillators and a placement of averaging operators

Research code and notes for the free oscillator

$$\ddot x+\omega_0^2x+\gamma|x|^\alpha\operatorname{sgn}x=0,\qquad x(0)=A,\ \dot x(0)=0,\qquad \alpha>0,\ \alpha\ne1 .$$

Every equivalent-linearization method that uses a fixed averaging operator gives
$\omega^2=\omega_0^2+\lambda\,\gamma A^{\alpha-1}$ with a **constant** $\lambda$. This includes the classical
(Krylov–Bogoliubov) average and the weighted averaging operators of Anh–Hai–Hieu (2017), Nguyen–Anh (2020) and
Nguyen et al. (2025). This repository asks where the **exact** value
$\lambda(q)=(\omega^2-\omega_0^2)/(\gamma A^{\alpha-1})$ lies, with $q=\gamma A^{\alpha-1}/\omega_0^2$ (the ratio of
the nonlinear to the linear restoring force at the amplitude), and places each operator relative to it.

## Main results

- **Exact representation.** With $x=A\sin\theta$, $1/\omega=\mathbb E_\theta[(1+qH(\theta))^{-1/2}]$, $\theta$ uniform on
  $[0,\pi/2]$, $H(\theta)=2(1-\sin^{\alpha+1}\theta)/((\alpha+1)\cos^2\theta)$ (a standard substitution).
- **Theorem (concavity and monotonicity).** $\omega^2$ is strictly concave in $q$ (Cauchy–Schwarz), so $\lambda(q)$
  decreases strictly from $\lambda_0=2\langle|\cos|^{\alpha+1}\rangle$ (the classical constant) to $\lambda_\infty$
  (the pure power-law oscillator), for every $\alpha>0$, $\alpha\ne1$, including $\alpha<1$ and softening springs.
  For the Duffing oscillator: $0.75\to0.7178$, initial slope $-3/128$.
- **Theorem (mixed power laws).**
  $\omega_0^2+\sum_i\lambda_{\infty,i}\gamma_iA^{\alpha_i-1}<\omega^2<\omega_0^2+\sum_i\lambda_{0,i}\gamma_iA^{\alpha_i-1}$.
- **Sign rule.** An operator with $\lambda_{\rm op}\ge\lambda_0$ always overestimates the frequency, one with
  $\lambda_{\rm op}\le\lambda_\infty$ always underestimates it, and one in between is exact at a single amplitude $q^\ast$.
- **Placement of the operators** for 14 values of $\alpha$ in $[0.1,11]$ (`output.txt`, `ranking.csv`): each operator is
  most accurate in its own range of $\alpha$, consistent with the examples chosen in the original papers.
- **Choice of weighting parameters** (`parameter_selection.py`): the minimax constant lies at 14–17 % of the
  bracket $(\lambda_\infty,\lambda_0)$ for every $\alpha$; the parameters that reproduce it are $s^\ast\approx2.05$
  (Anh 2017) and $p^\ast\approx0.48$ (Nguyen–Anh 2020) for the cubic case, close to the values used in those papers,
  and they vary with $\alpha$. The exponent $n$ of the superellipse family in Nguyen et al. (2025) barely changes $\lambda$.
- **Two-point Padé frequency formula** (`closed_form.py`): exact expansions at both ends ($q\to0$ to order $q^2$,
  $q\to\infty$ to order $1/q$) fix a [2/2] rational $\lambda(q)$ with no fitted parameter. Max frequency error:
  0.000048 % for Duffing (two AGM steps: 0.0014 %; rational harmonic balance, Beléndez et al. 2009: 0.0055 %),
  $\le0.012$ % for $\alpha\le11$, 0.054 % for $\alpha=20$. Five constants per $\alpha$ are tabulated in
  `closed_form_coefficients.csv`.
- **Power-mean linearization** (`power_mean_linearization.py`): keep the harmonic trial solution and the $x^2$ weight
  of classical linearization, but replace the arithmetic mean of the secant stiffness $g(x)/x$ by a power mean of
  order $r$. The optimal order is nearly universal, $r^\ast\in[0.437,0.442]$ for $0.1\le\alpha\le11$, and a single
  $r=0.44$ beats the best constant (minimax) linear operator for all 14 values of $\alpha$ (Duffing: 0.12 % vs 0.38 %).
- **Power-mean view** (`power_mean_view.py`): for any odd restoring force, with the energy chord stiffness $K(\theta)$,
  the exact $\omega^2$ is the power mean of order $-1/2$ of $K$, classical linearization is the arithmetic mean, and
  $M_{-1}(K)<\omega^2<M_0(K)<M_1(K)$. The parameter-free combination $\sqrt{M_0M_{-1}}$ is accurate to
  0.0014–0.09 % (Duffing, $x+x^7$, $x+x^3+x^5$, pendulum up to 3 rad, Duffing–harmonic); for Duffing it matches the
  accuracy of two AGM steps. Outside Duffing these means need one quadrature, as does the exact frequency. The value is conceptual: no linear averaging operator can be exact at all amplitudes.

The upper bound $\omega\le\sqrt{\omega_0^2+\lambda_0\gamma A^{\alpha-1}}$ (classical linearization / first-order
harmonic balance overestimates the frequency) is due to Benguria and Depassier (1996), for any odd nonlinearity.
The theorem above sharpens it: their bound is the $q\to0$ end of a monotone function whose other end is
$\lambda_\infty$. Rational harmonic balance (Beléndez et al. 2009) and the arithmetic–geometric mean give accurate
Duffing frequencies; the two-point Padé formula listed above improves on both and covers every $\alpha$.

Full statements, proofs, the complete table and the literature check are in [`NOTES_vi.md`](NOTES_vi.md)
(Vietnamese).

## Files

| File | Content |
|---|---|
| `frequency_bounds.py` | Checks the representation (elliptic integral, ODE), the limits and slope, monotonicity on a grid, the mixed and softening cases; places all operators; reproduces Tables 3 and 5 of Nguyen et al. (2025). About 15 s. |
| `output.txt`, `ranking.csv` | Output of the script above. |
| `NOTES_vi.md` | Notes with proofs, the full placement table and the literature check (Vietnamese). |
| `parameter_selection.py`, `parameter_selection.txt` | Minimax constant and optimal operator parameters $s^\ast$, $p^\ast$, $n^\ast$ for 14 values of $\alpha$. |
| `closed_form.py`, `closed_form.txt`, `closed_form_coefficients.csv` | Two-point Padé formula and its constants for 19 values of $\alpha$. |
| `power_mean_linearization.py`, `power_mean_linearization.txt` | Power-mean linearization: end-point orders $r_0$, $r_\infty$, optimal $r^\ast$ and the universal $r=0.44$. |
| `power_mean_view.py`, `power_mean_view.txt` | Power-mean bounds and parameter-free formulas for Duffing, $x+x^7$, $x+x^3+x^5$, pendulum, Duffing–harmonic. |
| `validation/duffing_weighted_averaging.py` | Cubic Duffing with the Anh–Hai–Hieu (2017) operator, $s=2$, against the exact frequency. |
| `validation/duffing_wao_2025.py` | Implementation of the Nguyen–Anh (2020) and Nguyen et al. (2025) operators; reproduces Table 7 of the 2020 paper and Tables 4–5 of the 2025 paper. |

## Running

```bash
pip install -r requirements.txt
python frequency_bounds.py
cd validation && python duffing_wao_2025.py
```

## References

1. N.D. Anh, N.Q. Hai, D.V. Hieu, The equivalent linearization method with a weighted averaging for analyzing of
   nonlinear vibrating systems, *Lat. Am. J. Solids Struct.* 14(9) (2017) 1723–1740. doi:10.1590/1679-78253488
2. Anh Tay Nguyen, N.D. Anh, A modified averaging operator with some applications, *Vietnam J. Mech.* 42(3) (2020)
   341–354. doi:10.15625/0866-7136/15356
3. Anh Tay Nguyen, N.N. Linh, N.C. Thang, N.A. Ngoc, L.Q. Vinh, N.D. Anh, A novel weighted averaging operator for
   frequency analysis of nonlinear free vibrations, *Int. J. Non-Linear Mech.* 173 (2025) 105057.
   doi:10.1016/j.ijnonlinmec.2025.105057
4. R.D. Benguria, M.C. Depassier, A variational principle for eigenvalue problems of Hamiltonian systems,
   *Phys. Rev. Lett.* 77 (1996) 2847–2850. doi:10.1103/PhysRevLett.77.2847
5. R.D. Benguria, M.C. Depassier, Variational calculation of the period of nonlinear oscillators, *J. Stat. Phys.*
   116 (2004) 923; arXiv:nlin/0409003.
6. J.D. García-Saldaña, A. Gasull, The period function and the harmonic balance method, *Bull. Sci. Math.* (2015);
   arXiv:1310.1215.
7. A. Beléndez, E. Gimeno, T. Beléndez, A. Hernández, Rational harmonic balance based method for conservative
   nonlinear oscillators: application to the Duffing equation, *Mech. Res. Commun.* 36(6) (2009) 728–734.
   doi:10.1016/j.mechrescom.2009.03.001

## Status

Research notes; a manuscript is in preparation.
