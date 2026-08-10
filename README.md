# THIẾT KẾ, MÔ PHỎNG VÀ ĐÓNG GÓI MÔ HÌNH BỘ LỌC THÔNG DẢI LC KHẢ CHỈNH BĂNG THÔNG KHÔNG ĐỔI DỰA TRÊN GHÉP HỖN HỢP ƯU THẾ TỪ TRÊN QUCS-S

> Phân tích lý thuyết → chứng minh mạng hai cổng → ma trận **Z** → ma trận **S** → kiểm chứng MATLAB/Qucs-S.

## 1. Mục tiêu và cơ sở tài liệu

Dự án xây dựng lại mô hình bộ lọc thông dải LC khả chỉnh dựa trên **magnetic-dominant mixed coupling**, sau đó chứng minh mô hình dưới dạng mạng hai cổng để thu được ma trận trở kháng **Z** và chuyển đổi sang ma trận tán xạ **S** dùng trong mô phỏng RF.

Tài liệu cơ sở chính:

> Longchuan Liu, Qian-Yin Xiang, Dinghong Jia, Xiaoguo Huang, and Quanyuan Feng,  
> **“A Novel Tunable LC Bandpass Filter with Constant Bandwidth based on Magnetic Dominant Mixed Coupling,”**  
> *Progress In Electromagnetics Research Letters*, Vol. 103, pp. 73–79, 2022.

Bài báo xác định:
- \(C_1\): phần tử ghép điện.
- \(L_4\): phần tử ghép từ.
- \(C_2\): tụ khả chỉnh dùng để điều chỉnh tần số trung tâm.
- Mạch được phân tích bằng hai chế độ **odd mode** và **even mode**.
- Bộ lọc thực nghiệm sử dụng varactor **SMV1212-004**.
- Kết quả đo của bài báo đạt dải điều chỉnh tần số trung tâm 72–222 MHz với băng thông \(-3\,\mathrm{dB}\) khoảng \(16.5\pm3.5\,\mathrm{MHz}\).

Phần **ma trận Z → ma trận S** trong repository này là phần mở rộng phân tích mạng hai cổng của dự án; bài báo gốc không trình bày trực tiếp phép chứng minh này.

---

## 2. Mạch, thông số và quy ước

Thông số dùng theo Figure 2 của bài báo:

| Tham số | Giá trị |
|---|---:|
| \(L_1\) | 82 nH |
| \(L_2\) | 82 nH |
| \(L_3\) | 68 nH |
| \(L_4\) | 18 nH |
| \(C_1\) | 0.2 pF |
| \(C_2\) | tụ varactor khả chỉnh |
| \(Z_0\) | 50 \(\Omega\) |

Trong phân tích miền Laplace:

$$
s=j\omega=j2\pi f
$$

Để biến đổi sang S-parameter chuẩn, hai dòng cổng phải được quy ước theo chuẩn mạng hai cổng: **dòng dương hướng vào mạng**. Nếu trên hình vẽ dòng \(I_{out}\) được chọn hướng ra khỏi cổng 2 thì cần đổi biến:

$$
I_2=-I_{out}
$$

trước khi đối chiếu trực tiếp với định nghĩa chuẩn của S-parameter.

Từ hệ phương trình vòng trên bảng, đặt:

$$
A=sL_3,\qquad B=sL_4,\qquad X=\frac{1}{sC_2}
$$

$$
P=s(L_1+L_3)
$$

$$
Q=s(L_2+L_3+L_4)+X
$$

$$
R=\frac{1}{sC_1}+2X
$$

Ba dòng vòng bên trong thỏa:

$$
\begin{bmatrix}
Q & -X & -B\\
-X & R & -X\\
-B & -X & Q
\end{bmatrix}
\begin{bmatrix}
i_2\\
i_3\\
i_4
\end{bmatrix}
=
\begin{bmatrix}
A i_1\\
0\\
-A I_{out}
\end{bmatrix}
$$

và hai điện áp cổng:

$$
U_1=P i_1-Ai_2
$$

$$
U_2=P I_{out}+Ai_4
$$

---

## 3. Chứng minh ma trận trở kháng Z và đối chiếu odd/even mode

Đặt:

$$
\mathbf{M}=
\begin{bmatrix}
Q & -X & -B\\
-X & R & -X\\
-B & -X & Q
\end{bmatrix}
$$

Định thức:

$$
\Delta=\det(\mathbf M)
=(Q+B)\left[R(Q-B)-2X^2\right]
$$

Hai phần tử nghịch đảo cần thiết:

$$
(\mathbf M^{-1})_{11}
=
\frac{QR-X^2}{\Delta}
$$

$$
(\mathbf M^{-1})_{13}
=
\frac{X^2+BR}{\Delta}
$$

Suy ra:

$$
Z_d
=
P-
A^2\frac{QR-X^2}{\Delta}
$$

$$
Z_m
=
A^2\frac{X^2+BR}{\Delta}
$$

Do cấu trúc đối xứng:

$$
\boxed{
\mathbf Z=
\begin{bmatrix}
Z_d & Z_m\\
Z_m & Z_d
\end{bmatrix}
}
$$

hay:

$$
Z_{11}=Z_{22}=Z_d,\qquad Z_{12}=Z_{21}=Z_m
$$

### Kiểm chứng bằng odd/even mode của bài báo

Bài báo cho admittance ngõ vào của odd mode và even mode:

$$
Y_{ino}
=
\frac{1}{
j\omega L_1+
\left[
j\omega L_3
\parallel
\left(
j\omega L_2+\frac{1}{j\omega(C_2+2C_1)}
\right)
\right]
}
$$

$$
Y_{ine}
=
\frac{1}{
j\omega L_1+
\left[
j\omega L_3
\parallel
\left(
j\omega(L_2+2L_4)+\frac{1}{j\omega C_2}
\right)
\right]
}
$$

Do đó:

$$
Z_o=\frac{1}{Y_{ino}},\qquad Z_e=\frac{1}{Y_{ine}}
$$

Với mạng hai cổng đối xứng:

$$
Z_e=Z_{11}+Z_{12}
$$

$$
Z_o=Z_{11}-Z_{12}
$$

nên:

$$
\boxed{
Z_{11}=Z_{22}=\frac{Z_e+Z_o}{2}
}
$$

$$
\boxed{
Z_{12}=Z_{21}=\frac{Z_e-Z_o}{2}
}
$$

Đây là phép kiểm chứng độc lập cho kết quả thu được từ hệ phương trình vòng trên bảng.

---

## 4. Chuyển đổi ma trận Z sang ma trận S

Với hai cổng cùng trở kháng tham chiếu:

$$
Z_{01}=Z_{02}=Z_0=50\,\Omega
$$

công thức ma trận là:

$$
\boxed{
\mathbf S=
(\mathbf Z-Z_0\mathbf I)
(\mathbf Z+Z_0\mathbf I)^{-1}
}
$$

Đặt:

$$
\Delta_S=
(Z_{11}+50)(Z_{22}+50)-Z_{12}Z_{21}
$$

ta được:

$$
S_{11}
=
\frac{
(Z_{11}-50)(Z_{22}+50)-Z_{12}Z_{21}
}{
\Delta_S
}
$$

$$
S_{12}
=
\frac{100Z_{12}}{\Delta_S}
$$

$$
S_{21}
=
\frac{100Z_{21}}{\Delta_S}
$$

$$
S_{22}
=
\frac{
(Z_{11}+50)(Z_{22}-50)-Z_{12}Z_{21}
}{
\Delta_S
}
$$

Với mạch đối xứng và thuận nghịch:

$$
S_{11}=S_{22},\qquad S_{12}=S_{21}
$$

### Hình ma trận S dùng trực tiếp trong GitHub

![Ma trận S từ ma trận Z](assets/s_matrix_github.png)

Hình trên dùng kiểu chữ toán học Computer Modern, phù hợp với cách GitHub/KaTeX hiển thị biểu thức toán.

---

## 5. Kiểm chứng tại 100 MHz và hướng triển khai Qucs-S

Từ điều kiện cộng hưởng của \(Y_{ino}\) và \(Y_{ine}\), có thể rút ra:

$$
f_o=
\frac{1}{2\pi}
\sqrt{
\frac{L_1+L_3}
{
(C_2+2C_1)
\left[
L_1(L_2+L_3)+L_2L_3
\right]
}
}
$$

$$
f_e=
\frac{1}{2\pi}
\sqrt{
\frac{L_1+L_3}
{
C_2
\left[
L_1(L_2+L_3+2L_4)+L_3(L_2+2L_4)
\right]
}
}
$$

và:

$$
f_c=\frac{f_o+f_e}{2}
$$

Với các giá trị linh kiện của Figure 2 và yêu cầu:

$$
f_c=100\,\mathrm{MHz}
$$

nghiệm của mô hình LC lý tưởng là xấp xỉ:

$$
\boxed{C_2\approx18.50\,\mathrm{pF}}
$$

Khi đó:

$$
f_o\approx106.06\,\mathrm{MHz}
$$

$$
f_e\approx93.94\,\mathrm{MHz}
$$

và tại \(100\,\mathrm{MHz}\):

$$
\mathbf Z
\approx
\begin{bmatrix}
-j26.423 & j61.444\\
j61.444 & -j26.423
\end{bmatrix}\Omega
$$

Suy ra:

$$
\mathbf S
\approx
\begin{bmatrix}
0.08452+j0.04004 & -0.42628+j0.89974\\
-0.42628+j0.89974 & 0.08452+j0.04004
\end{bmatrix}
$$

tương ứng:

$$
|S_{11}|\approx-20.58\,\mathrm{dB}
$$

$$
|S_{21}|\approx-0.038\,\mathrm{dB}
$$

Đây là **kết quả của mô hình LC lý tưởng**, chưa bao gồm tổn hao cuộn cảm, ESR/ký sinh của varactor, tụ ghép chéo \(C_3\), đường truyền PCB và tổn hao điện môi.

Hướng triển khai trên Qucs-S:

`Mạch LC lý tưởng → sweep C2 → Z-parameter → S-parameter → varactor thực → ký sinh linh kiện → so sánh với bài báo → đóng gói thành subcircuit/model reusable`

### Phân biệt phần kế thừa và phần thực hiện của dự án

**Kế thừa từ bài báo**
- Topology magnetic-dominant mixed coupling.
- Phân tích odd/even mode.
- Giá trị \(L_1,L_2,L_3,L_4,C_1\).
- Varactor SMV1212-004 và nguyên lý điều chỉnh \(C_2\).
- Tiêu chí đánh giá theo S-parameter và băng thông.

**Phần triển khai của dự án**
- Thiết lập hệ KVL từ mạch.
- Chứng minh ma trận \(Z\).
- Chuyển đổi \(Z\rightarrow S\).
- Sweep \(C_2\) tại các tần số mục tiêu.
- Mô phỏng và đóng gói mô hình trên Qucs-S.
- Đối chiếu mô hình lý tưởng với mô hình linh kiện thực.

---

### Reference

L. Liu, Q.-Y. Xiang, D. Jia, X. Huang, and Q. Feng, “A Novel Tunable LC Bandpass Filter with Constant Bandwidth based on Magnetic Dominant Mixed Coupling,” *Progress In Electromagnetics Research Letters*, vol. 103, pp. 73–79, 2022.
