# THIẾT KẾ, MÔ PHỎNG VÀ ĐÓNG GÓI MÔ HÌNH BỘ LỌC THÔNG DẢI LC KHẢ CHỈNH BĂNG THÔNG KHÔNG ĐỔI DỰA TRÊN GHÉP HỖN HỢP ƯU THẾ TỪ TRÊN QUCS-S

**Mục tiêu:** xây dựng lại mạch LC khả chỉnh từ tài liệu tham khảo, chứng minh mô hình mạng hai cổng bằng hệ phương trình mạch, suy ra ma trận trở kháng **Z**, chuyển đổi sang ma trận tán xạ **S**, sau đó dùng kết quả làm cơ sở kiểm chứng trên **MATLAB** và **Qucs-S**.

> **Lưu ý về hiển thị GitHub:** toàn bộ công thức nhiều dòng trong README này dùng khối `math` của GitHub thay cho việc đặt `$$` ngay sau tiêu đề. Cách này tránh lỗi xuống dòng của `bmatrix`, dấu `_`, và dấu `\\` trong ma trận.

---

## 1. Cơ sở tài liệu và thông số mạch

Tài liệu tham khảo chính:

**Longchuan Liu, Qian-Yin Xiang, Dinghong Jia, Xiaoguo Huang, and Quanyuan Feng**,  
**“A Novel Tunable LC Bandpass Filter with Constant Bandwidth based on Magnetic Dominant Mixed Coupling,”**  
*Progress In Electromagnetics Research Letters*, Vol. 103, pp. 73–79, 2022.

Bài báo sử dụng cấu trúc **magnetic-dominant mixed coupling**. Trong cấu trúc này:

- $C_1$ là phần tử ghép điện.
- $L_4$ là phần tử ghép từ.
- $C_2$ là tụ khả chỉnh dùng để điều chỉnh tần số trung tâm.
- Hai mạch tương đương **odd mode** và **even mode** được dùng để phân tích cộng hưởng và hệ số ghép.
- Mạch thực nghiệm của bài báo sử dụng varactor **SMV1212-004**.

Các giá trị dùng trong phần phân tích của dự án:

| Tham số | Giá trị |
|---|---:|
| $L_1$ | 82 nH |
| $L_2$ | 82 nH |
| $L_3$ | 68 nH |
| $L_4$ | 18 nH |
| $C_1$ | 0.2 pF |
| $C_2$ | Varactor, biến thiên |
| $Z_0$ | 50 Ω |

Trong miền tần số:

```math
s=j\omega=j2\pi f
```

Phần **Z → S** trong repository này là phần chứng minh mở rộng của dự án. Bài báo gốc tập trung vào odd/even mode, coupling coefficient, mô phỏng và đo kiểm; bài báo không trình bày trực tiếp phép suy ra ma trận **Z** từ hệ KVL như phần dưới đây.

---

## 2. Chứng minh từ hệ phương trình trên bảng

Ảnh dưới đây là hệ phương trình mạch được dùng làm điểm xuất phát:

![Hệ phương trình chứng minh trên bảng](assets/whiteboard_derivation.png)

Đặt các đại lượng trung gian:

```math
A=sL_3
```

```math
B=sL_4
```

```math
X=\frac{1}{sC_2}
```

```math
P=s(L_1+L_3)
```

```math
Q=s(L_2+L_3+L_4)+X
```

```math
R=\frac{1}{sC_1}+2X
```

Từ ba phương trình dòng điện nội, hệ được viết gọn thành:

```math
\begin{bmatrix}
Q & -X & -B \\
-X & R & -X \\
-B & -X & Q
\end{bmatrix}
\begin{bmatrix}
i_2 \\
i_3 \\
i_4
\end{bmatrix}
=
\begin{bmatrix}
A i_1 \\
0 \\
-A I_{\mathrm{out}}
\end{bmatrix}
```

Đặt:

```math
\mathbf{M}=
\begin{bmatrix}
Q & -X & -B \\
-X & R & -X \\
-B & -X & Q
\end{bmatrix}
```

Khi đó:

```math
\begin{bmatrix}
i_2 \\
i_3 \\
i_4
\end{bmatrix}
=
\mathbf{M}^{-1}
\begin{bmatrix}
A i_1 \\
0 \\
-A I_{\mathrm{out}}
\end{bmatrix}
```

Định thức của $\mathbf M$:

```math
\Delta
=
\det(\mathbf{M})
=
(Q+B)\left[R(Q-B)-2X^2\right]
```

Hai phần tử cần dùng của $\mathbf M^{-1}$:

```math
\left(\mathbf{M}^{-1}\right)_{11}
=
\frac{QR-X^2}{\Delta}
```

```math
\left(\mathbf{M}^{-1}\right)_{13}
=
\frac{X^2+BR}{\Delta}
```

Từ đó suy ra hai thành phần trở kháng đặc trưng của mạng đối xứng:

```math
Z_d
=
P
-
A^2
\frac{QR-X^2}{\Delta}
```

```math
Z_m
=
A^2
\frac{X^2+BR}{\Delta}
```

Ma trận trở kháng có dạng:

```math
\boxed{
\mathbf{Z}
=
\begin{bmatrix}
Z_d & Z_m \\
Z_m & Z_d
\end{bmatrix}
}
```

hay:

```math
Z_{11}=Z_{22}=Z_d
```

```math
Z_{12}=Z_{21}=Z_m
```

### Quy ước dòng cổng

Khi đưa kết quả sang định nghĩa S-parameter chuẩn, dòng ở cả hai port phải được quy ước **dương khi đi vào mạng hai cổng**. Nếu mũi tên $I_{\mathrm{out}}$ trên sơ đồ được vẽ hướng ra khỏi port 2, cần đổi dấu dòng port 2 trước khi lập ma trận Z chuẩn.

---

## 3. Kiểm chứng bằng odd mode và even mode của bài báo

Bài báo cho admittance ở odd mode:

```math
Y_{\mathrm{ino}}
=
\frac{1}{
j\omega L_1
+
\left[
j\omega L_3
\parallel
\left(
j\omega L_2
+
\frac{1}{j\omega(C_2+2C_1)}
\right)
\right]
}
```

và admittance ở even mode:

```math
Y_{\mathrm{ine}}
=
\frac{1}{
j\omega L_1
+
\left[
j\omega L_3
\parallel
\left(
j\omega(L_2+2L_4)
+
\frac{1}{j\omega C_2}
\right)
\right]
}
```

Suy ra:

```math
Z_o=\frac{1}{Y_{\mathrm{ino}}}
```

```math
Z_e=\frac{1}{Y_{\mathrm{ine}}}
```

Với một mạng hai cổng đối xứng:

```math
Z_e=Z_{11}+Z_{12}
```

```math
Z_o=Z_{11}-Z_{12}
```

Do đó:

```math
\boxed{
Z_{11}=Z_{22}=\frac{Z_e+Z_o}{2}
}
```

```math
\boxed{
Z_{12}=Z_{21}=\frac{Z_e-Z_o}{2}
}
```

Kết quả này là bước kiểm chứng độc lập cho ma trận **Z** suy ra từ hệ phương trình trên bảng.

Tần số cộng hưởng của hai mode được xác định từ điều kiện của bài báo:

```math
\operatorname{Im}\left(Y_{\mathrm{ino}}\right)=0
```

```math
\operatorname{Im}\left(Y_{\mathrm{ine}}\right)=0
```

Sau khi rút gọn, tần số odd mode:

```math
f_o
=
\frac{1}{2\pi}
\sqrt{
\frac{L_1+L_3}{
(C_2+2C_1)
\left[
L_1(L_2+L_3)+L_2L_3
\right]
}
}
```

Tần số even mode:

```math
f_e
=
\frac{1}{2\pi}
\sqrt{
\frac{L_1+L_3}{
C_2
\left[
L_1(L_2+L_3+2L_4)
+
L_3(L_2+2L_4)
\right]
}
}
```

Tần số trung tâm:

```math
f_c=\frac{f_o+f_e}{2}
```

Với bộ giá trị linh kiện ở trên và yêu cầu:

```math
f_c=100\ \mathrm{MHz}
```

mô hình LC lý tưởng cho nghiệm xấp xỉ:

```math
\boxed{
C_2\approx18.50\ \mathrm{pF}
}
```

tương ứng:

```math
f_o\approx106.06\ \mathrm{MHz}
```

```math
f_e\approx93.94\ \mathrm{MHz}
```

---

## 4. Chuyển đổi ma trận Z sang ma trận S

Hai cổng được chuẩn hóa cùng trở kháng:

```math
Z_{01}=Z_{02}=Z_0=50\ \Omega
```

Công thức ma trận chuẩn:

```math
\boxed{
\mathbf{S}
=
\left(
\mathbf{Z}-Z_0\mathbf{I}
\right)
\left(
\mathbf{Z}+Z_0\mathbf{I}
\right)^{-1}
}
```

Với $Z_0=50\ \Omega$:

```math
\boxed{
\mathbf{S}
=
\left(
\mathbf{Z}-50\mathbf{I}
\right)
\left(
\mathbf{Z}+50\mathbf{I}
\right)^{-1}
}
```

Đặt:

```math
\Delta_S
=
(Z_{11}+50)(Z_{22}+50)
-
Z_{12}Z_{21}
```

Các phần tử S-parameter:

```math
S_{11}
=
\frac{
(Z_{11}-50)(Z_{22}+50)
-
Z_{12}Z_{21}
}{
\Delta_S
}
```

```math
S_{12}
=
\frac{
100Z_{12}
}{
\Delta_S
}
```

```math
S_{21}
=
\frac{
100Z_{21}
}{
\Delta_S
}
```

```math
S_{22}
=
\frac{
(Z_{11}+50)(Z_{22}-50)
-
Z_{12}Z_{21}
}{
\Delta_S
}
```

Do mạng đối xứng và thuận nghịch:

```math
S_{11}=S_{22}
```

```math
S_{12}=S_{21}
```

Ma trận S có thể viết:

```math
\boxed{
\mathbf{S}
=
\begin{bmatrix}
S_{11} & S_{12} \\
S_{21} & S_{22}
\end{bmatrix}
}
```

### Hình ma trận S dùng trong GitHub

![Quan hệ ma trận Z sang ma trận S](assets/s_matrix_github.png)

Nếu GitHub hiển thị SVG tốt hơn trên repository, có thể thay đường dẫn ảnh bằng:

```text
assets/s_matrix_github.svg
```

### Kiểm tra số tại 100 MHz

Với $C_2\approx18.4966\ \mathrm{pF}$:

```math
\mathbf{Z}
\approx
\begin{bmatrix}
-j26.423 & j61.444 \\
j61.444 & -j26.423
\end{bmatrix}
\ \Omega
```

Suy ra:

```math
\mathbf{S}
\approx
\begin{bmatrix}
0.08452+j0.04004 & -0.42628+j0.89974 \\
-0.42628+j0.89974 & 0.08452+j0.04004
\end{bmatrix}
```

Về biên độ:

```math
20\log_{10}|S_{11}|
\approx
-20.58\ \mathrm{dB}
```

```math
20\log_{10}|S_{21}|
\approx
-0.038\ \mathrm{dB}
```

Các giá trị này thuộc **mô hình LC lý tưởng** tại 100 MHz, chưa phải kết quả của mạch thực nghiệm hoàn chỉnh trong bài báo.

---

## 5. Hướng triển khai và đóng gói mô hình trên Qucs-S

Chuỗi triển khai đề xuất:

```text
Mạch LC từ bài báo
        ↓
Kiểm tra odd mode / even mode
        ↓
Thiết lập hệ phương trình mạch
        ↓
Suy ra ma trận Z
        ↓
Chuyển Z → S với Z0 = 50 Ω
        ↓
Sweep C2 và tần số
        ↓
So sánh MATLAB ↔ Qucs-S
        ↓
Thay C2 lý tưởng bằng mô hình varactor
        ↓
Bổ sung tổn hao và ký sinh linh kiện
        ↓
Đóng gói thành subcircuit/model tái sử dụng trong Qucs-S
```

### Phần kế thừa từ bài báo

- Topology bộ lọc LC khả chỉnh.
- Cấu trúc **magnetic-dominant mixed coupling**.
- Odd-mode và even-mode equivalent circuits.
- Các giá trị $L_1$, $L_2$, $L_3$, $L_4$, $C_1$ dùng trong Figure 2.
- Varactor SMV1212-004 trong mạch thực nghiệm.
- Tiêu chí đánh giá theo $S_{11}$, $S_{21}$, tần số trung tâm và băng thông.

### Phần thực hiện trong dự án

- Viết lại hệ phương trình mạch từ topology.
- Chứng minh ma trận **Z**.
- Kiểm chứng ma trận Z bằng odd/even mode.
- Chuyển đổi **Z → S**.
- Xây dựng chương trình MATLAB để sweep $C_2$.
- Dựng mạch và kiểm chứng trên Qucs-S.
- Bổ sung mô hình linh kiện thực và ký sinh.
- Đóng gói mạch thành subcircuit/model dùng lại.

### Giới hạn của mô hình hiện tại

Mô hình lý tưởng chưa bao gồm:

- ESR và hệ số $Q$ hữu hạn của cuộn cảm.
- Điện trở nối tiếp và ký sinh của varactor.
- Tụ ghép chéo $C_3$ của mạch hoàn chỉnh trong bài báo.
- Đường truyền PCB và discontinuity.
- Tổn hao điện môi và tổn hao dẫn.
- Sai số linh kiện.

Vì vậy, kết quả MATLAB của mô hình LC lý tưởng dùng để **xác minh phần chứng minh**, không được đồng nhất trực tiếp với kết quả đo của prototype trong bài báo.

---

### Tài liệu tham khảo chính

L. Liu, Q.-Y. Xiang, D. Jia, X. Huang, and Q. Feng,  
“A Novel Tunable LC Bandpass Filter with Constant Bandwidth based on Magnetic Dominant Mixed Coupling,”  
*Progress In Electromagnetics Research Letters*, Vol. 103, pp. 73–79, 2022.
