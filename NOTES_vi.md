# Chặn hai phía cho tần số dao động tự do và xếp chỗ các toán tử trung bình

Ghi chú nghiên cứu (tiếng Việt). Ý chính: mỗi phương pháp tuyến tính hóa với một toán tử trung bình cố định
là một hằng số; đặt hằng số đó so với khoảng chặn của hằng số đúng để biết dấu và độ lớn sai số.
Bài toán: tần số dao động tự do theo biên độ.

- Script: `frequency_bounds.py` (chạy khoảng 15 s). Kết quả: `output.txt`, `ranking.csv`.
- Chứng minh ở **mục 1–3**, xếp chỗ các toán tử ở **mục 4–5**, kiểm tra tài liệu ở **mục 6**.

## 0. Bài toán và ký hiệu

$$
\ddot x+\omega_0^2x+\gamma|x|^\alpha\operatorname{sgn}x=0,\qquad x(0)=A,\ \dot x(0)=0,\qquad
\alpha>0,\ \alpha\ne1 .
$$
Lấy $\omega_0=1$ (đổi thang thời gian), đặt $q=\gamma A^{\alpha-1}/\omega_0^2$ và viết
$$
\omega^2=\omega_0^2+\lambda(q)\,\gamma A^{\alpha-1}\qquad\Big(\lambda(q)=\frac{\omega^2-1}{q}\Big).
$$
Mọi phương pháp tuyến tính hóa với một toán tử trung bình cố định (cổ điển, Anh–Hai–Hieu 2017,
Nguyen–Anh 2020, Nguyen et al. 2025) đều cho **một hằng số** $\lambda_{\rm op}$, không phụ thuộc $A,\gamma$,
vì tỉ số $\langle|\cos|^{\alpha+1}\rangle_{\rm op}/\langle\cos^2\rangle_{\rm op}$ không chứa biên độ.
Vì vậy câu hỏi "toán tử nào tốt hơn" quy về câu hỏi: hằng số đúng $\lambda(q)$ nằm ở đâu và thay đổi
thế nào theo $q$.

## 1. Biểu diễn chính xác (Bổ đề 1)

Đặt
$$
h(s)=\frac{2\,(1-s^{\alpha+1})}{(\alpha+1)(1-s^2)},\quad h(1)=1,\qquad
H(\theta)=h(\sin\theta),\qquad \Theta\sim\text{đều trên }[0,\tfrac\pi2].
$$
**Bổ đề 1.** Với mọi $q>-1/\max H$ (lò xo cứng, hoặc lò xo mềm trước khi mất dao động),
$$
\frac1\omega=\mathbb E\big[(1+qH)^{-1/2}\big],\qquad\text{tức là}\qquad
\omega^2=F(q):=\Big(\mathbb E\big[(1+qH)^{-1/2}\big]\Big)^{-2}.
$$

*Chứng minh.* Bảo toàn năng lượng cho $T=4\int_0^A dx/\sqrt{2(U(A)-U(x))}$, với
$U(x)=x^2/2+\gamma|x|^{\alpha+1}/(\alpha+1)$. Đổi biến $x=As$:
$$
U(A)-U(As)=A^2\Big[\frac{1-s^2}2+\frac{q(1-s^{\alpha+1})}{\alpha+1}\Big]=\frac{A^2(1-s^2)}2\big[1+q\,h(s)\big].
$$
Đặt tiếp $s=\sin\theta$, ta có $ds/\sqrt{1-s^2}=d\theta$, nên
$T=4\int_0^{\pi/2}(1+qH)^{-1/2}d\theta=2\pi\,\mathbb E[(1+qH)^{-1/2}]$ và $\omega=2\pi/T$. $\square$

Hàm $h$ liên tục và dương trên $[0,1]$. Nó **không hằng** khi $\alpha\ne1$, vì $h(0)=2/(\alpha+1)\ne1=h(1)$.

*Kiểm số:* khớp công thức elliptic ($\alpha=3$) tới $10^{-11}$, khớp tích phân ODE ($\alpha=0{,}5;\,2{,}5;\,7$)
tới $10^{-12}$.

## 2. Định lý: $\lambda(q)$ giảm ngặt (Định lý 2)

**Định lý 2.** Với mọi $\alpha>0$, $\alpha\ne1$, hàm $F(q)=\omega^2$ **lõm ngặt** trên $(-1/\max H,\infty)$. Do đó:

1. $\lambda(q)=\big(F(q)-F(0)\big)/q$ **giảm ngặt** trên toàn miền, kể cả $q<0$;
2. $\lambda(0^+)=\lambda_0:=\mathbb E[H]=2\langle|\cos\theta|^{\alpha+1}\rangle$, đúng bằng hằng số của
   **tuyến tính hóa cổ điển** (Krylov–Bogoliubov, cũng là cân bằng điều hòa bậc nhất);
3. $\lambda(\infty)=\lambda_\infty:=\big(\mathbb E[H^{-1/2}]\big)^{-2}$, là bình phương tần số của dao động
   **thuần lũy thừa** $\ddot x+|x|^\alpha\operatorname{sgn}x=0$ với $A=1$;
4. $\lambda'(0)=-\tfrac34\operatorname{Var}(H)$; với Duffing, $\lambda'(0)=-3/128$, khớp Lindstedt.

Như vậy $\lambda_\infty<\lambda(q)<\lambda_0$ với mọi $q>0$.

*Chứng minh.* Đặt $Y=1+qH>0$ và $m(q)=\mathbb E[Y^{-1/2}]$. Khi đó $m'=-\tfrac12\mathbb E[HY^{-3/2}]$ và
$m''=\tfrac34\mathbb E[H^2Y^{-5/2}]$. Với $F=m^{-2}$:
$$
F''=2m^{-4}\big(3m'^2-m\,m''\big)
=\tfrac32m^{-4}\Big(\big(\mathbb E[HY^{-3/2}]\big)^2-\mathbb E[Y^{-1/2}]\,\mathbb E[H^2Y^{-5/2}]\Big)\le0
$$
theo bất đẳng thức Cauchy–Schwarz cho $Y^{-1/4}$ và $HY^{-5/4}$. Dấu bằng chỉ xảy ra khi $H\propto Y=1+qH$,
tức $H$ là hằng, điều đã bị loại. Vậy $F$ lõm ngặt, và độ dốc dây cung $\lambda(q)=(F(q)-F(0))/q$
của một hàm lõm ngặt thì giảm ngặt.

Về các giới hạn:
- $F'(0)=-2m'(0)=\mathbb E H$.
- $F(q)/q=\big(\mathbb E[(1/q+H)^{-1/2}]\big)^{-2}\to(\mathbb E H^{-1/2})^{-2}$ theo hội tụ trội (vì $H\ge\min h>0$),
  nên $\lambda=F/q-1/q\to\lambda_\infty$.
- Khai triển $F=1+\mathbb E[H]\,q-\tfrac34\operatorname{Var}(H)\,q^2+O(q^3)$ cho độ dốc ban đầu.

Đẳng thức $\mathbb E H=2\langle|\cos|^{\alpha+1}\rangle$ có được bằng tích phân từng phần với $d(\tan\theta)$:
$$
\int_0^{\pi/2}\frac{1-\sin^{\alpha+1}\theta}{\cos^2\theta}d\theta=(\alpha+1)\int_0^{\pi/2}\sin^{\alpha+1}\theta\,d\theta ,
$$
vì số hạng biên bằng $0$ ở cả hai đầu. Cuối cùng, $\lambda_\infty$ trùng tần số dao động thuần lũy thừa: đặt $q\to\infty$
trong Bổ đề 1, hoặc so trực tiếp với công thức Gamma
$\omega=\sqrt{\pi(\alpha+1)/2}\,\Gamma\big(\tfrac{\alpha+3}{2(\alpha+1)}\big)/\Gamma\big(\tfrac1{\alpha+1}\big)$
(khớp tới $10^{-10}$). $\square$

**Hệ quả 1 (quy tắc dấu).** Với hằng số $\lambda_c$, sai số
$\omega_c^2-\omega^2=q\,(\lambda_c-\lambda(q))$ có dấu xác định hoàn toàn bởi vị trí của $\lambda_c$:
- $\lambda_c\ge\lambda_0$: luôn đánh giá **cao** tần số;
- $\lambda_c\le\lambda_\infty$: luôn đánh giá **thấp**;
- $\lambda_\infty<\lambda_c<\lambda_0$: đúng tại **đúng một** $q^\ast$, đánh giá cao khi $q<q^\ast$ và thấp khi $q>q^\ast$.

**Hệ quả 2 (cổ điển luôn đánh giá cao).** Tuyến tính hóa cổ điển đánh giá cao tần số **ở mọi biên độ**, với mọi
$\alpha>0$, $\alpha\ne1$. Điều này đúng cả cho lò xo cứng ($q>0$) lẫn lò xo mềm ($q<0$, trước khi mất dao động),
vì với $q<0$ thì $\lambda(q)>\lambda_0$ và $q\lambda(q)<q\lambda_0$. Kiểm số ở $\alpha=3$, $q=-0{,}5$:
sai số $+0{,}77\,\%$.

> **Không mới.** Hệ quả 2 đã có trong Benguria & Depassier (PRL 77, 1996, Eq. (18)–(19)), với mọi phi tuyến
> lẻ $N(u)$ và mọi biên độ. Họ dùng một nguyên lý biến phân với hàm thử tối ưu cho bài toán tuyến tính, và
> cận thu được chính là $\omega\le\sqrt{1+\tfrac34\delta a^2}$ cho Duffing. Khi viết bài phải trích dẫn họ cho
> kết quả này, và chỉ coi Hệ quả 2 là hệ quả phụ của Định lý 2.

*Kiểm số Định lý 2:* $\lambda(q)$ giảm ngặt trên 400 điểm $10^{-3}\le q\le10^6$, với 14 giá trị
$\alpha\in[0{,}1;\,11]$.

## 3. Định lý: chặn hai phía cho lực đa thức hỗn hợp (Định lý 3)

Xét $g(x)=\omega_0^2x+\sum_i\gamma_i|x|^{\alpha_i}\operatorname{sgn}x$, với $\gamma_i>0$, $\alpha_i>0$, $\alpha_i\ne1$,
và $q_i=\gamma_iA^{\alpha_i-1}/\omega_0^2$. Bổ đề 1 vẫn đúng với $Y=1+\sum_iq_iH_i$.

**Định lý 3.**
$$
\omega_0^2+\sum_i\lambda_{\infty,i}\,\gamma_iA^{\alpha_i-1}\;<\;\omega^2\;<\;\omega_0^2+\sum_i\lambda_{0,i}\,\gamma_iA^{\alpha_i-1}.
$$
Cận trên chính là tuyến tính hóa cổ điển. Cận dưới bằng $\omega_0^2$ cộng tổng các bình phương tần số của
từng dao động thuần lũy thừa.

*Chứng minh.* $M(Y)=(\mathbb E\,Y^{-1/2})^{-2}$ là trung bình lũy thừa bậc $-1/2$.
- Cận trên: $M(Y)<\mathbb E Y$ (bất đẳng thức trung bình lũy thừa, chặt vì $Y$ không hằng).
- Cận dưới: phép tính $F''$ ở Định lý 2, với $H$ thay bằng một hướng bất kỳ $Z$, cho thấy $M$ lõm.
  $M$ cũng thuần nhất bậc 1, nên **siêu cộng tính**: $M(X+Z)=2M\big(\tfrac{X+Z}2\big)\ge M(X)+M(Z)$, với dấu
  bằng chỉ khi $X\propto Z$. Áp dụng lặp, $M(Y)>M(1)+\sum_iM(q_iH_i)=1+\sum_iq_i\lambda_{\infty,i}$. $\square$

*Kiểm số:* $x+q_3x^3+q_5x^5$, với $0\le q_3\le100$ và $10^{-3}\le q_5\le10^3$: đúng ở 25/25 điểm.

## 4. Xếp chỗ khách quan bốn toán tử (phần c)

Các toán tử được cài đúng theo định nghĩa trong bài gốc:
- **cổ điển:** trung bình đều trên một chu kỳ;
- **Anh–Hai–Hieu 2017** (LAJSS 14:1723): $\int_0^\infty s^2\theta e^{-s\theta}(\cdot)d\theta$ với $s=2$;
- **Nguyen–Anh 2020** (VJMech 42:341, WLAO): $p$ cố định $0{,}25;\,0{,}5;\,1$, lấy min theo $h$;
- **Nguyen et al. 2025** (IJNLM 173:105057, WEL): $p=\pm\sqrt{1-h^2}$, lấy min theo $h$ và dấu.

Cài đặt được kiểm bằng cách tái lập Bảng 3 (EL và WEL, $\alpha=0{,}5;\,3;\,10$) và Bảng 5 (cột $\omega_{AW}$, $n=1..5$) của
bài 2025, cùng Bảng 7 của bài 2020 (xem thư mục `validation/`).

**Đọc bảng.** Mỗi ô ghi ký hiệu dấu và sai số tần số lớn nhất, tính bằng % trên $10^{-3}\le q\le10^6$:
- ↑: luôn đánh giá cao ($\lambda_{\rm op}\ge\lambda_0$);
- ↓: luôn đánh giá thấp ($\lambda_{\rm op}\le\lambda_\infty$);
- ×$q^\ast$: chính xác tại $q^\ast$, đổi dấu một lần.

Hai cột cuối không phải toán tử của bài nào. "$\lambda_\infty$" là cận dưới dùng như một hằng số. "Hữu tỷ"
là hằng số phụ thuộc mức (mục 5).

| α | λ₀ | λ∞ | Cổ điển | Anh 2017 s=2 | WLAO p=0,25 | WLAO p=0,5 | WLAO p=1 | WEL 2025 | λ∞ (cận) | Hữu tỷ |
|---|---|---|---|---|---|---|---|---|---|---|
| 0,1 | 1,2358 | 1,2062 | ↑ 1,22 | ↑ 1,82 | ×0,76 0,60 | ×37 0,28 | ↓ 1,03 | ↓ 1,02 | 0,29 | 0,00 |
| 0,25 | 1,1852 | 1,1668 | ↑ 0,79 | ↑ 1,32 | ×1,3 0,29 | ↓ 0,28 | ↓ 1,01 | ↓ 1,00 | 0,19 | 0,00 |
| 0,5 | 1,1128 | 1,1059 | ↑ 0,31 | ↑ 0,70 | ×470 0,08 | ↓ 0,29 | ↓ 0,81 | ↓ 0,81 | 0,08 | 0,00 |
| 0,75 | 1,0520 | 1,0505 | ↑ 0,07 | ↑ 0,28 | ↓ 0,08 | ↓ 0,21 | ↓ 0,46 | ↓ 0,46 | 0,02 | 0,00 |
| 1,5 | 0,9153 | 0,9114 | ↑ 0,21 | ↓ 0,25 | ↓ 0,09 | ↓ 0,37 | ↓ 1,05 | ↓ 1,04 | 0,05 | 0,00 |
| 2 | 0,8488 | 0,8366 | ↑ 0,73 | ↓ 0,33 | ×3 0,21 | ↓ 0,40 | ↓ 1,62 | ↓ 1,61 | 0,18 | 0,00 |
| 2,5 | 0,7949 | 0,7728 | ↑ 1,42 | ↓ 0,41 | ×1,4 0,69 | ↓ 0,41 | ↓ 1,89 | ↓ 1,88 | 0,36 | 0,00 |
| 3 | 0,7500 | 0,7178 | ↑ 2,22 | ×19 0,48 | ×0,99 1,30 | ×8,8 0,41 | ↓ 1,97 | ↓ 1,95 | 0,56 | 0,01 |
| 4 | 0,6791 | 0,6279 | ↑ 3,99 | ×6,3 0,80 | ×0,71 2,75 | ×3 1,39 | ↓ 2,02 | ↓ 2,00 | 0,99 | 0,02 |
| 5 | 0,6250 | 0,5578 | ↑ 5,86 | ×4,9 1,53 | ×0,59 4,32 | ×2,1 2,65 | ↓ 2,05 | ↓ 2,03 | 1,42 | 0,05 |
| 7 | 0,5469 | 0,4556 | ↑ 9,56 | ×4,2 3,02 | ×0,5 7,54 | ×1,5 5,35 | ×61 2,07 | ×55 2,05 | 2,21 | 0,11 |
| 9 | 0,4922 | 0,3849 | ↑ 13,08 | ×4,2 4,42 | ×0,45 10,66 | ×1,3 8,03 | ×13 2,06 | ×12 2,05 | 2,90 | 0,20 |
| 10 | 0,4703 | 0,3571 | ↑ 14,76 | ×4,3 5,08 | ×0,44 12,16 | ×1,2 9,33 | ×9,8 2,82 | ×9,6 2,85 | 3,20 | 0,25 |
| 11 | 0,4512 | 0,3331 | ↑ 16,38 | ×4,3 5,70 | ×0,43 13,61 | ×1,2 10,59 | ×8,3 3,66 | ×8,2 3,69 | 3,48 | 0,30 |

**Nhận xét (mô tả, không đánh giá):**
1. Mỗi toán tử có một vùng $\alpha$ mà ở đó nó cho sai số lớn nhất nhỏ nhất trong sáu cột toán tử:
   - WLAO $p=0{,}5$ khi $\alpha\le0{,}25$;
   - WLAO $p=0{,}25$ khi $\alpha=0{,}5;\,1{,}5;\,2$;
   - cổ điển và WLAO $p=0{,}25$ gần bằng nhau khi $\alpha=0{,}75$;
   - WLAO $p=0{,}5$ khi $\alpha=2{,}5$–$3$ (ngang Anh 2017 ở $2{,}5$);
   - Anh 2017 ($s=2$) khi $\alpha=4$–$5$;
   - WEL 2025 và WLAO $p=1$ khi $\alpha\ge7$.

   Điều này phù hợp với các ví dụ mà từng bài gốc đã chọn để minh họa.
2. Toán tử nào có $\lambda_{\rm op}$ nằm trong khoảng $(\lambda_\infty,\lambda_0)$ thì chính xác tại đúng một
   biên độ $q^\ast$. Vị trí $q^\ast$ cho biết toán tử đó được "căn chỉnh" cho vùng phi tuyến nào. Với $\alpha\ge3$:
   - WLAO $p=0{,}25$: $q^\ast\approx0{,}4$–$1$ (phi tuyến vừa);
   - Anh 2017: $q^\ast\approx4$–$19$;
   - WEL 2025 (khi $\alpha\ge7$): $q^\ast\approx8$–$55$ (phi tuyến mạnh).
3. WEL 2025 và WLAO 2020 với $p=1$ cho hằng số gần như trùng nhau (lệch dưới $0{,}1\,\%$) ở mọi $\alpha$. Lý do là
   cực tiểu của WEL nằm ở $h^\ast\approx0{,}07$–$0{,}09$, nơi $p=-\sqrt{1-h^2}\approx-1$.
4. Các toán tử có $\lambda_{\rm op}\le\lambda_\infty$ (WEL 2025 và WLAO $p=1$ khi $\alpha\le5$, cùng các ô ↓ khác)
   đánh giá thấp tần số ở mọi biên độ. Trong các hằng số có tính chất này, $\lambda_\infty$ có sai số nhỏ nhất
   ở mọi biên độ, theo Hệ quả 1.

## 5. Hằng số phụ thuộc mức qua hai giới hạn đúng

Nội suy hữu tỷ qua hai giới hạn đúng và độ dốc ban đầu:
$$
\lambda(q)\approx\frac{\lambda_0+\lambda_\infty b\,q}{1+b\,q},\qquad
b=\frac{\tfrac34\operatorname{Var}(H)}{\lambda_0-\lambda_\infty}.
$$
Công thức này khớp đúng ba thông tin đã chứng minh: $\lambda_0$, $\lambda_\infty$ và $\lambda'(0)$.
- Duffing: $b=0{,}727$, sai số tần số nằm trong $-0{,}0074\,\%$ đến $0$.
- Mọi $\alpha\le11$: sai số không quá $0{,}30\,\%$, và luôn **đánh giá thấp**.

Tính chất luôn đánh giá thấp mới chỉ là quan sát số, chưa chứng minh.

**Lưu ý quan trọng (xem mục 6):** với riêng Duffing, công thức này **không** phải chính xác nhất.
- Hai bước trung bình cộng–nhân (AGM) của Gauss cho
  $\omega\approx\sqrt{1+q}\,\big(1+(1-m)^{1/4}\big)^2/4$, với $m=q/(2(1+q))$. Sai số lớn nhất là $0{,}0014\,\%$
  (đã tự kiểm lại), không cần tham số khớp nào.
- Cân bằng điều hòa hữu tỷ bậc hai (Beléndez et al. 2009) báo sai số khoảng $0{,}0055\,\%$ khi $A\to\infty$.

Ưu điểm duy nhất còn lại của công thức hữu tỷ là nó dùng được cho **mọi** $\alpha$ (AGM chỉ áp dụng cho bậc 3).
Ưu điểm này chưa được so sánh với các phương pháp tổng quát khác, nên **không nên** dùng "độ chính xác cao"
làm điểm bán của bài.

## 6. Kiểm tra tài liệu (phần a)

Bước đầu tra cứu qua kết quả tìm kiếm web. Sau khi mở mạng, **toàn văn** bốn bài quan trọng nhất đã được đọc (PDF
không kèm trong repo vì bản quyền). Ba bài lấy từ arXiv; bài Beléndez 2009 lấy bản preprint
từ kho RUA của Đại học Alicante, vì bản xuất bản bị khóa. Các bài còn lại trong bảng mới chỉ xem tóm tắt.

### 6.1 Công thức hữu tỷ hai đầu (mục 5): đã có kết quả tốt hơn

| Tài liệu | Nội dung | Quan hệ |
|---|---|---|
| Beléndez, Gimeno, Beléndez, Hernández, *Mech. Res. Commun.* 36 (2009) 728–734, doi 10.1016/j.mechrescom.2009.03.001 (**đã đọc**, Bảng 1) | cân bằng điều hòa hữu tỷ bậc hai cho Duffing: sai số $0{,}0015\,\%$ ($\lambda=5$), $0{,}0027\,\%$ ($10$), $0{,}0051\,\%$ ($100$), $0{,}0055\,\%$ ($\infty$); xét cả lò xo mềm | chính xác hơn công thức của mình ($0{,}0074\,\%$) |
| Beléndez et al., *J. Appl. Math.* 2012, 286290, doi 10.1155/2012/286290; Carvalhaes & Suppes, *Am. J. Phys.* 76 (2008) 1150 | tần số qua AGM bằng hàm sơ cấp | công thức AGM 2 bước: $0{,}0014\,\%$ (đã tự kiểm) |
| Beléndez et al., *Phys. Lett. A* (2009), phương pháp lập phương hóa; Elías-Zúñiga et al., *Math. Probl. Eng.* 2012, 618750 | công thức khớp cả hai giới hạn, sai số 0,07 % và 0,055 % | tương tự ý tưởng hai đầu |
| Lai, Lim, Wu et al., *Appl. Math. Model.* (2009) | Newton–cân bằng điều hòa cho Duffing bậc 3–5 | tham khảo |
| Cveticanin, Zuković, Cveticanin, *Mathematics* 12 (2024) 3040 | tần số dựng từ tần số chính xác của dao động thuần bậc ba | gần với ý dùng $\lambda_\infty$ |

**Kết luận:** chưa thấy bài nào dùng **đúng** dạng hữu tỷ với $\lambda_0$, $\lambda_\infty$ và độ dốc chính xác,
nhưng về độ chính xác nó thua các kết quả đã có. Mục 5 chỉ nên giữ như một minh họa.

### 6.2 Định lý 2–3: một phần đã có, phần chính có vẻ mới

| Tài liệu (đã đọc toàn văn) | Nội dung thực tế | Quan hệ với kết quả của mình |
|---|---|---|
| Benguria & Depassier, *Phys. Rev. Lett.* 77 (1996) 2847 | Nguyên lý biến phân $\lambda=\max_g(\dots)$ cho bài toán trị riêng phi tuyến. Với hàm thử $g=u/\sqrt{1-u^2/u_m^2}$, Eq. (18) cho cận áp dụng **với mọi** $N(u)$; với Duffing, Eq. (19) cho $\omega\le\sqrt{1+\tfrac34\delta a^2}$ ở mọi biên độ. | **Trùng Hệ quả 2** (cổ điển hay cân bằng điều hòa luôn đánh giá cao tần số), và họ làm tổng quát hơn. Bài không có tính lõm, không có $\lambda(q)$ đơn điệu, không có cận dưới. |
| Benguria & Depassier, *J. Stat. Phys.* 116 (2004) 923 | Nguyên lý đối ngẫu (Fenchel–Moreau) cho **cận dưới** của $\omega$, Eq. (31)–(32). Hàm thử là khai triển Taylor của $(1+\delta z)^{-1/2}$, nên các cận dưới trùng với các đường nhiễu loạn (tác giả tự nhận xét như vậy). | Có ý "chặn hai phía", nhưng cận dưới của họ là phân thức đa thức có hệ số bậc cao nhất âm, nên **mất ý nghĩa khi biên độ lớn**. Cận dưới $\lambda_\infty$ của mình đơn giản, sắc ở biên độ lớn và đúng ở mọi biên độ: hai bên bổ sung cho nhau. |
| García-Saldaña & Gasull, *Bull. Sci. Math.* (2015) | Chứng minh hàm chu kỳ $T(A)$ đơn điệu cho vài họ thế; cân bằng điều hòa bậc 1–3 giữ được tính đơn điệu, dáng điệu tại $0$ và tại $\infty$. Có dùng phép đổi biến $z=x/A$ giống Bổ đề 1. | Không có định lý về dấu sai số, tính lõm hay $\lambda$ đơn điệu. Bổ đề 1 là phép đổi biến chuẩn, không mới. |
| Opial (1961); Chicone (1987); Schaaf (1985); Chow & Sanders (1986); Cima, Gasull, Mañosas (2000); Benguria, Depassier, Loss (2016) (chỉ xem tóm tắt) | tính đơn điệu của **hàm chu kỳ** $T(A)$ | yếu hơn: $\omega(A)$ đơn điệu không suy ra $\lambda(q)$ đơn điệu |
| Gottlieb, *J. Sound Vib.* 261 (2003) 557 (chỉ xem tóm tắt) | so cân bằng điều hòa với nghiệm đúng cho lũy thừa phân số | tham khảo |

**Kết luận sau khi đọc toàn văn:**
- **Đã có:**
  - Bổ đề 1 (phép đổi biến $x=A\sin\theta$, chuẩn trong tài liệu);
  - Hệ quả 2 (Benguria & Depassier 1996, tổng quát hơn);
  - ý tưởng có cận dưới (Benguria & Depassier 2004, nhưng là cận kiểu nhiễu loạn).
- **Chưa thấy ở đâu:**
  1. $\omega^2$ **lõm ngặt** theo $q$ và $\lambda(q)$ **giảm ngặt** từ $\lambda_0$ xuống $\lambda_\infty$, với mọi $\alpha>0$, $\alpha\ne1$ (Định lý 2), chứng minh gọn bằng Cauchy–Schwarz;
  2. **cận dưới** $\omega^2>\omega_0^2+\lambda_\infty\gamma A^{\alpha-1}$, dạng đóng, sắc khi $A\to\infty$, và bản siêu cộng tính cho đa thức hỗn hợp (Định lý 3);
  3. **quy tắc dấu** (Hệ quả 1) và **bảng xếp chỗ** các toán tử trung bình hằng số của nhóm GS. Anh (mục 4).
- **Cách viết đề xuất:** trình bày Định lý 2 như sự làm sắc kết quả của Benguria & Depassier 1996, cụ thể là
  "cận trên của họ là đầu $q\to0$ của một hàm đơn điệu, mà đầu kia là $\lambda_\infty$". Trọng tâm bài đặt vào
  quy tắc dấu và việc xếp chỗ các toán tử.

### 6.3 Toán tử trung bình (để trích dẫn đúng)

- Anh, Hai, Hieu, *Lat. Am. J. Solids Struct.* 14 (2017) 1723–1740.
- Hieu, Hai, Hung, *J. Appl. Math.* 2018, 7487851, doi 10.1155/2018/7487851 (tác giả là Hieu, Hai, **Hung**).
- Anh Tay Nguyen, N.D. Anh, *Vietnam J. Mech.* 42(3) (2020) 341–354, doi 10.15625/0866-7136/15356.
- Anh Tay Nguyen, N.N. Linh, N.C. Thang, N.A. Ngoc, L.Q. Vinh, N.D. Anh, *Int. J. Non-Linear Mech.* 173 (2025)
  105057, doi 10.1016/j.ijnonlinmec.2025.105057.
- N.D. Anh, *Vietnam J. Mech.* (2015), cách nhìn đối ngẫu về giá trị trung bình, dạng của hệ số trọng số.

## 7. Giới hạn và việc còn mở

- Định lý 2 và 3 cần lực phục hồi dạng lũy thừa có tâm, không cản. Hệ có cản hoặc lực phụ thuộc vận tốc
  nằm ngoài phạm vi.
- Hằng số của WLAO và WEL được tìm bằng cách quét $h$ (4001 điểm, sau đó Brent). Chữ số thứ ba có thể lệch
  so với bài gốc, vì bài gốc quét $h$ thô hơn.
- Tính luôn-đánh-giá-thấp của hằng số hữu tỷ ở mục 5 mới là quan sát số.
