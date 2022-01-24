#include<iostream>
#include<cstdio>
#include<algorithm>
#include<cmath>
#include<vector>
#include<cstring> 
using namespace std;

//定义复数
struct Complex {
	double Vr, Vi;
	Complex(double Ur = 0, double Ui = 0) :Vr(Ur), Vi(Ui) {}
};

//复数的四则运算
Complex operator +(Complex Va, Complex Vb) { return Complex(Va.Vr + Vb.Vr, Va.Vi + Vb.Vi); }
Complex operator -(Complex Va, Complex Vb) { return Complex(Va.Vr - Vb.Vr, Va.Vi - Vb.Vi); }
Complex operator *(Complex Va, Complex Vb) { return Complex(Va.Vr * Vb.Vr - Va.Vi * Vb.Vi, Va.Vr * Vb.Vi + Va.Vi * Vb.Vr); }
Complex operator /(Complex Va, double Vb) { return Complex(Va.Vr / Vb, Va.Vi / Vb); }

typedef Complex cplx;
const double PI = acos(-1);   //π
const int N = 1e6+10;
const int t = 7;     //明文模数
const int q = 874;   //密文模数（系数模）
const int q1 = -q / 2 + 1;
const int q2 = q / 2;
const int n = 16;    //向量维度
int k = q / t;

vector<int> s, a, e, pk0, pk1;   //私钥，用于生成公钥的多项式a，噪音多项式e，公钥
vector<int> ct0, ct1;   //将密文定义为全局变量，便于传参
cplx Ca[N], Cb[N], Cc[N];



//对多项式的系数和次数进行取模，系数模q，次数模n，A的低位对应多项式的低次幂
vector<int> Mod(vector<int> A) {
	int L = A.size();
	//次数取模
	for (int i = L - 1; i >= n; i--) {   //i为次数
		A[i - n] = A[i - n] - A[i];
	}
	//系数取模
	for (int i = 0; i < n; i++) {
		while (A[i] > q2 )   A[i] = A[i] - q;
		while (A[i] < q1 )   A[i] = A[i] + q;
	}
	vector<int> B(n);
	for (int i = 0; i < n; i++)   B[i] = A[i];
	return B;
}



//求多项式的点值表示，参数Ci为多项式的系数表示，执行Change函数之后Ci为点值表示
void Change(cplx* Ci, int len, int tmp) {
	int lglen = log2(len);
	for (int i = 0; i < len; i++) {
		int Revi = 0;
		for (int j = 0; j < lglen; j++) Revi |= (((i >> j) & 1) << (lglen - 1 - j));
		if (Revi < i) swap(Ci[Revi], Ci[i]);
	}
	for (int s = 2; s <= len; s <<= 1) {
		int t = s / 2;
		cplx Vper{ cos(tmp * PI / t),sin(tmp * PI / t) };
		for (int i = 0; i < len; i += s) {
			cplx Vr{ 1,0 };
			for (int j = 0; j < t; j++, Vr = Vr * Vper) {
				cplx Ve = Ci[i + j], Vo = Ci[i + j + t];
				Ci[i + j] = Ve + Vo * Vr;
				Ci[i + j + t] = Ve - Vo * Vr;
			}
		}
	}
	if (tmp == -1)
		for (int i = 0; i <= len; i++)
			Ci[i] = Ci[i] / (double)len;
}



vector<int> FFT(vector<int> A, vector<int> B) {
    memset(Ca,0,sizeof(struct Complex)*N);
    memset(Cb,0,sizeof(struct Complex)*N);
    memset(Cc,0,sizeof(struct Complex)*N);
    
	for (int i = 0; i < n; i++)   Ca[i].Vr = A[i];   //多项式系数
	for (int i = 0; i < n; i++)   Cb[i].Vr = B[i];

	int Nc = 2 * (n - 1);   //乘积的幂次
	int len = 1 << (int)ceil(log2(Nc));   //ceil:取整函数
	if (len == Nc) len <<= 1;
	Change(Ca, len, 1);
	Change(Cb, len, 1);
	for (int i = 0; i <= len; i++) Cc[i] = Ca[i] * Cb[i];
	Change(Cc, len, -1);   //逆FFT

	//保存结果到mul
	vector<int> mul;
	for (int i = 0; i <= Nc; i++)   mul.push_back((int)round(Cc[i].Vr));
	
	return Mod(mul);
}



//不同长度的两个多项式相加
vector<int> Add(vector<int> A, vector<int> B) {
	int L = max(A.size(), B.size());
	vector<int> C(L, 0);
	for (int i = 0; i < L; i++) {
		if (i >= A.size())   C[i] = B[i];
		else if (i >= B.size())   C[i] = A[i];
		else   C[i] = A[i] + B[i];
	}
	return Mod(C);
}



//随机生成私钥和公钥
void init() {
	//随机生成私钥s
	for (int i = n - 1; i >= 0; i--) {
		int sig1 = rand();   //±1
		if (sig1 % 2)   s.push_back(rand() % 2);
		else s.push_back(-(rand() % 2));
	}
	printf("随机生成的私钥s为：\n");
	for (int i = 0; i < s.size(); i++) {
		printf("%d ", s[i]);
		if (i == s.size() - 1)   printf("\n");
	}

	//随机生成a，也就是pk1
	for (int i = n - 1; i >= 0; i--) {
		int sig2 = rand();
		if (sig2 % 2)   a.push_back(rand() % q2);
		else   a.push_back(-rand() % q2);
	}
	pk1 = a;
	
	//生成噪音e
	for (int i = n - 1; i >= 0; i--) {
		int sig3 = rand();
		if (sig3 % 2)   e.push_back(rand() % 10);
		else   e.push_back(-rand() % 10);
	}

	//计算a*s
	vector<int> as = FFT(a, s);
	//计算-as
	int L = as.size();
	for (int i = 0; i < L; i++)   as[i] = -as[i];
	//计算-as+e
	pk0 = Add(as, e);
}



void Enc(vector<int> m) {
	//随机生成e1,e2,u
	vector<int> e1, e2, u;
	for (int i = n - 1; i >= 0; i--) {
		int sig4 = rand();
		if (sig4 % 2)   e1.push_back(rand() % 10);
		else   e1.push_back(-rand() % 10);
	}
	for (int i = n - 1; i >= 0; i--) {
		int sig5 = rand();
		if (sig5 % 2)   e2.push_back(rand() % 10);
		else   e2.push_back(-rand() % 10);
	}
	for (int i = n - 1; i >= 0; i--) {
		int sig6 = rand();   //±1
		if (sig6 % 2)  u.push_back(rand() % 2);
		else u.push_back(-(rand() % 2));
	}

	//计算密文ct
	for (int i = 0; i < m.size(); i++)   m[i] *= k;   //k=q/t

	vector<int> pk0u = FFT(pk0, u);
	ct0 = Add(Add(pk0u, e1), m);
	vector<int> pk1u = FFT(pk1, u);
	ct1 = Add(pk1u, e2);
	printf("对明文m加密后，得到密文ct0和ct1\n");
	printf("ct0为：");
	for (int i = 0; i < ct0.size(); i++)   printf("%d ", ct0[i]);
	printf("\nct1为：");
	for (int i = 0; i < ct1.size(); i++)   printf("%d ", ct1[i]);
}



void Dec(vector<int> ct0, vector<int> ct1) {   //传入密文参数ct
	vector<int> M = Add(ct0, FFT(ct1, s));
	printf("\n解密出的明文M为:\n");
	for (int i = 0; i < M.size(); i++) {
		M[i] /= k;
		printf("%d ", M[i]);
	}
}



int main() {
	init();
	vector<int> m = {3,2,1};
	printf("明文消息m为:\n");
	for (int i = m.size() - 1; i >= 0; i--) {
		printf("%d ", m[i]);
		if (i == 0)   printf("\n");
	}
	Enc(m);
	Dec(ct0, ct1);
	return 0;
}
