#ifndef LOW_LEVEL_FFT_H
#define LOW_LEVEL_FFT_H

/*!
Fast Fourier/Cosine/Sine Transform	<BR>
dimension   :one <BR>
data length :power of 2<BR>
decimation  :frequency<BR>
radix       :split-radix<BR>
data        :inplace<BR>
table       :not use<BR>
functions<BR>
cdft: Complex Discrete Fourier Transform<BR>
rdft: Real Discrete Fourier Transform<BR>
ddct: Discrete Cosine Transform<BR>
ddst: Discrete Sine Transform<BR>
dfct: Cosine Transform of RDFT (Real Symmetric DFT)<BR>
dfst: Sine Transform of RDFT (Real Anti-symmetric DFT)<BR>
macro definitions<BR>
USE_CDFT_PTHREADS : default=not defined<BR>
CDFT_THREADS_BEGIN_N  : must be >= 512, default=8192<BR>
CDFT_4THREADS_BEGIN_N : must be >= 512, default=65536<BR>
USE_CDFT_WINTHREADS : default=not defined<BR>
CDFT_THREADS_BEGIN_N  : must be >= 512, default=32768<BR>
CDFT_4THREADS_BEGIN_N : must be >= 512, default=524288<BR>


-------- Complex DFT (Discrete Fourier Transform) --------<BR>
X[k] = sum_j=0^n-1 x[j]*exp(2*pi*i*j*k/n), 0<=k<n<BR>
X[k] = sum_j=0^n-1 x[j]*exp(-2*pi*i*j*k/n), 0<=k<n<BR>
(notes: sum_j=0^n-1 is a summation from j=0 to n-1)<BR>
cdft(2*n, 1, a);<BR>
cdft(2*n, -1, a);<BR>
2*n            :data length (int)<BR>
n >= 1, n = power of 2
a[0...2*n-1]   :input/output data (double *)<BR>
input data<BR>
a[2*j] = Re(x[j]), <BR>
a[2*j+1] = Im(x[j]), 0<=j<n<BR>
output data<BR>
a[2*k] = Re(X[k]), <BR>
a[2*k+1] = Im(X[k]), 0<=k<n<BR>
Inverse of <BR>
cdft(2*n, -1, a);<BR>
is <BR>
cdft(2*n, 1, a);<BR>
for (j = 0; j <= 2 * n - 1; j++) {<BR>
a[j] *= 1.0 / n;<BR>
<BR>

-------- Real DFT / Inverse of Real DFT --------<BR>
[definition]
<case1> RDFT
R[k] = sum_j=0^n-1 a[j]*cos(2*pi*j*k/n), 0<=k<=n/2
I[k] = sum_j=0^n-1 a[j]*sin(2*pi*j*k/n), 0<k<n/2
<case2> IRDFT (excluding scale)
a[k] = (R[0] + R[n/2]*cos(pi*k))/2 + 
sum_j=1^n/2-1 R[j]*cos(2*pi*j*k/n) + 
sum_j=1^n/2-1 I[j]*sin(2*pi*j*k/n), 0<=k<n
[usage]
<case1>
rdft(n, 1, a);
<case2>
rdft(n, -1, a);
[parameters]
n              :data length (int)
n >= 2, n = power of 2
a[0...n-1]     :input/output data (double *)
<case1>
output data
a[2*k] = R[k], 0<=k<n/2
a[2*k+1] = I[k], 0<k<n/2
a[1] = R[n/2]
<case2>
input data
a[2*j] = R[j], 0<=j<n/2
a[2*j+1] = I[j], 0<j<n/2
a[1] = R[n/2]
[remark]
Inverse of 
rdft(n, 1, a);
is 
rdft(n, -1, a);
for (j = 0; j <= n - 1; j++) {
a[j] *= 2.0 / n;
}
.


-------- DCT (Discrete Cosine Transform) / Inverse of DCT --------
[definition]
<case1> IDCT (excluding scale)
C[k] = sum_j=0^n-1 a[j]*cos(pi*j*(k+1/2)/n), 0<=k<n
<case2> DCT
C[k] = sum_j=0^n-1 a[j]*cos(pi*(j+1/2)*k/n), 0<=k<n
[usage]
<case1>
ddct(n, 1, a);
<case2>
ddct(n, -1, a);
[parameters]
n              :data length (int)
n >= 2, n = power of 2
a[0...n-1]     :input/output data (double *)
output data
a[k] = C[k], 0<=k<n
[remark]
Inverse of 
ddct(n, -1, a);
is 
a[0] *= 0.5;
ddct(n, 1, a);
for (j = 0; j <= n - 1; j++) {
a[j] *= 2.0 / n;
}
.


-------- DST (Discrete Sine Transform) / Inverse of DST --------
[definition]
<case1> IDST (excluding scale)
S[k] = sum_j=1^n A[j]*sin(pi*j*(k+1/2)/n), 0<=k<n
<case2> DST
S[k] = sum_j=0^n-1 a[j]*sin(pi*(j+1/2)*k/n), 0<k<=n
[usage]
<case1>
ddst(n, 1, a);
<case2>
ddst(n, -1, a);
[parameters]
n              :data length (int)
n >= 2, n = power of 2
a[0...n-1]     :input/output data (double *)
<case1>
input data
a[j] = A[j], 0<j<n
a[0] = A[n]
output data
a[k] = S[k], 0<=k<n
<case2>
output data
a[k] = S[k], 0<k<n
a[0] = S[n]
[remark]
Inverse of 
ddst(n, -1, a);
is 
a[0] *= 0.5;
ddst(n, 1, a);
for (j = 0; j <= n - 1; j++) {
a[j] *= 2.0 / n;
}
.


-------- Cosine Transform of RDFT (Real Symmetric DFT) --------
[definition]
C[k] = sum_j=0^n a[j]*cos(pi*j*k/n), 0<=k<=n
[usage]
dfct(n, a);
[parameters]
n              :data length - 1 (int)
n >= 2, n = power of 2
a[0...n]       :input/output data (double *)
output data
a[k] = C[k], 0<=k<=n
[remark]
Inverse of 
a[0] *= 0.5;
a[n] *= 0.5;
dfct(n, a);
is 
a[0] *= 0.5;
a[n] *= 0.5;
dfct(n, a);
for (j = 0; j <= n; j++) {
a[j] *= 2.0 / n;
}
.


-------- Sine Transform of RDFT (Real Anti-symmetric DFT) --------
[definition]
S[k] = sum_j=1^n-1 a[j]*sin(pi*j*k/n), 0<k<n
[usage]
dfst(n, a);
[parameters]
n              :data length + 1 (int)
n >= 2, n = power of 2
a[0...n-1]     :input/output data (double *)
output data
a[k] = S[k], 0<k<n
(a[0] is used for work area)
[remark]
Inverse of 
dfst(n, a);
is 
dfst(n, a);
for (j = 1; j <= n - 1; j++) {
a[j] *= 2.0 / n;
}
.
*/
class LowLevelFFT
{
public:
	LowLevelFFT();
	~LowLevelFFT();
public:
	static void cdft(int, int, double *);
	static void rdft(int, int, double *);
	static void ddct(int, int, double *);
	static void ddst(int, int, double *);
	static void dfct(int, double *);
	static void dfst(int, double *);
private:
	inline static void cftfsub(int n, double *a);
	inline static void cftbsub(int n, double *a);
	inline static void bitrv2(int n, double *a);
	inline static void bitrv2conj(int n, double *a);
	inline static void bitrv216(double *a);
	inline static void bitrv216neg(double *a);
	inline static void bitrv208(double *a);
	inline static void bitrv208neg(double *a);
	inline static void bitrv1(int n, double *a);
	inline static void cftb1st(int n, double *a);
	inline static void cftrec4_th(int n, double *a);
	inline static void *cftrec1_th(void *p);
	inline static void *cftrec2_th(void *p);
	inline static void cftrec4(int n, double *a);
	inline static int cfttree(int n, int j, int k, double *a);
	inline static void cftleaf(int n, int isplt, double *a);
	inline static void cftmdl1(int n, double *a);
	inline static void cftmdl2(int n, double *a);
	inline static void cftfx41(int n, double *a);
	inline static void cftf161(double *a);
	inline static void cftf162(double *a);
	inline static void cftf081(double *a);
	inline static void cftf082(double *a);
	inline static void cftf040(double *a);
	inline static void cftb040(double *a);
	inline static void cftx020(double *a);
	inline static void rftfsub(int n, double *a);
	inline static void rftbsub(int n, double *a);
	inline static void dctsub(int n, double *a);
	inline static void dstsub(int n, double *a);
	inline static void dctsub4(int n, double *a);
	inline static void dstsub4(int n, double *a);
};

#endif
