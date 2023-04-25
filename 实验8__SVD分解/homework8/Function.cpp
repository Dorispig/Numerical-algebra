#include<iostream>
#include <vector>
#include"Function.h"
#include<math.h>
using namespace std;
vector<vector<double>> Matrix;
//矩阵相加
vector<vector<double>> Matrix_add(vector<vector<double>> a, vector<vector<double>> b) {
	for (int i = 0; i < a.size(); i++) {
		for (int j = 0; j < a[0].size(); j++)
			a[i][j] += b[i][j];
	}
	return a;
}
//矩阵相乘
vector<vector<double>> Matrix_Mul(vector<vector<double>> a, vector<vector<double>> b) {
	vector<vector<double>> p;
	int m, n, s;
	m = a.size();
	n = a[0].size();
	s = b[0].size();
	p = make_zero_matrix(m,s);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < s; j++) {
			for (int k = 0; k < n; k++) 
				p[i][j] += a[i][k] * b[k][j];
		}
	}
	return p;
}
//矩阵乘向量
vector<double> Matrix_Mul_Vector(vector<vector<double>> A, vector<double> x) {
	double s;
	int m, n;
	m = A.size(); n = A[0].size();
	vector<double> b(m, 0);
	for (int i = 0; i < m; i++) {
		s = 0;
		for (int j = 0; j < n; j++)
			s += A[i][j] * x[j];
		b[i] = s;
	}
	return b;
}
//向量乘矩阵
vector<double> Vector_Mul_Matrix(vector<double> x, vector<vector<double>> A) {
	int m, n;
	m = A.size(); n = A[0].size();
	vector<double> b(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			b[i] += x[j] * A[j][i];
	}
	return b;
}
//列向量乘行向量
vector<vector<double>> LVector_Mul_HVector(vector<double> x, vector<double> y) {
	vector<vector<double>> A;
	int m, n;
	m = x.size(); n = y.size();
	A = make_zero_matrix(m, n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			A[i][j] = x[i] * y[j];
	}
	return A;
}
//行向量乘列向量
double HVector_Mul_LVector(vector<double> x, vector<double> y) {
	int n = x.size();
	double t = 0;
	for (int i = 0; i < n; i++)
		t += x[i]*y[i];
	return t;
}
//求矩阵转置
vector<vector<double>> Matrix_T(vector<vector<double>> a) {
	vector<vector<double>> t;
	int m, n;
	m = a.size(); n = a[0].size();
	t = make_zero_matrix(n,m);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			t[j][i] = a[i][j];
	}
	return t;
}
//矩阵数乘
vector<vector<double>> Matrix_Mul_num(vector<vector<double>> a, double b) {
	for (int i = 0; i < a.size(); i++) {
		for (int j = 0; j < a[0].size(); j++)
			a[i][j] *= b;
	}
	return a;
}
//向量数乘
vector<double> Vector_Mul_num(vector<double> a, double b) {
	for (int j = 0; j < a.size(); j++)
		a[j] *= b;
	return a;
}
//矩阵除以一个数
vector<vector<double>> Matrix_Divid_num(vector<vector<double>> a, double b) {
	for (int i = 0; i < a.size(); i++) {
		for (int j = 0; j < a[0].size(); j++)
			a[i][j] /= b;
	}
	return a;
}
//矩阵相减,a-b
vector<vector<double>> Matrix_Subtracition(vector<vector<double>> a, vector<vector<double>> b) {
	for (int i = 0; i < a.size(); i++) {
		for (int j = 0; j < a[0].size(); j++)
			a[i][j] -= b[i][j];
	}
	return a;
}
//向量相减
vector<double> Vector_Subtracition(vector<double> a, vector<double> b) {
	vector<double> c(a.size());
	for (int i = 0; i < a.size(); i++)
			c[i] = a[i] - b[i];
	return c;
} 
//向量相加
vector<double> Vector_Add(vector<double> a, vector<double> b) {
	vector<double> c(a.size());
	for (int i = 0; i < a.size(); i++)
		c[i] = a[i] + b[i];
	return c;
}
//创建0矩阵
vector<vector<double>> make_zero_matrix(int a, int b) {
	vector<double> m(b);
	vector<vector<double>> n(a, m);
	return n;
}
//化为上三角矩阵求det（A）
double det(vector<vector<double>> A) {
	double m=1;
	int n = A.size();
	if (n == 1)
		return A[0][0];
	else {
			for (int i = 0; i < n - 1; i++) {
				for (int j = i + 1; j < n; j++) {
					for (int k = i; k < n; k++)
						A[j][k] -= A[i][k] * A[j][i] / A[i][i];
				}
			}
			for (int i = 0; i < n; i++)
				m *= A[i][i];
			return m;
	}	
}
//不选主元的Gauss消去法求L,U,直接用for循环
vector<double> Gauss_no(vector<vector<double>> A,vector<double> b){
	vector<vector<double>>L;
	vector<vector<double>>U;
	int n = A.size();
	L = make_zero_matrix(n, n);
	U = make_zero_matrix(n, n);
	for (int k = 0; k < n-1; k++) {
		for (int i = k + 1; i < n; i++) {
			A[i][k] /= A[k][k];
			for (int j = k + 1; j < n; j++)
				A[i][j] -= A[i][k] * A[k][j];
		}
	}
	getLU(A, L, U);
	b = Lower_triangle( L, b);
	b = Upper_triangle( U, b);
	return b;
}
//取L,U矩阵
void getLU(vector<vector<double>> A, vector<vector<double>>& L, vector<vector<double>>& U) {
	int n = A.size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				L[i][j] = 1;
				U[i][j] = A[i][j];
			}
			else if (i < j) U[i][j] = A[i][j];
			else L[i][j] = A[i][j];
		}
	}
}
//分块矩阵(从start_x,start_y开始，大小为block_x * block_y的分块矩阵)
vector<vector<double>> Matrix_block(vector<vector<double>> A, int start_x, int start_y, int block_x, int block_y) {
	vector<vector<double>> B;
	B = make_zero_matrix(block_x, block_y);
	for (int i = 0; i < block_x; i++) {
		for (int j = 0; j < block_y; j++) {
			B[i][j] = A[i+ start_x][j + start_y];
		}
	}
	return B;
}
//将A中分块矩阵更新为乘B之后的值,k=-1左乘B,k=1右乘B
vector<vector<double>> Matrix_block_mul(vector<vector<double>> A, int start_x, int start_y, int block_x, int block_y, vector<vector<double>>B,int k) {
	vector<vector<double>> C= Matrix_block(A, start_x, start_y, block_x, block_y),D;
	D = make_zero_matrix(block_x, block_y);
	D = Matrix_block(A, start_x, start_y, block_x, block_y);
	if (k == 1) D = Matrix_Mul(D, B); 
	else D = Matrix_Mul(B, D);
	for (int i = start_x; i < start_x + block_x; i++) {
		for (int j = start_y; j < start_y + block_y; j++) A[i][j] = D[i - start_x][j - start_y];
	}
	return A;
}
//将A中分块矩阵更新为加B之后的值,k=-1减k=1加
vector<vector<double>> Matrix_block_add(vector<vector<double>> A, int start_x, int start_y, int block_x, int block_y, vector<vector<double>>B, int k) {
	for (int i = start_x; i < start_x + block_x; i++) {
		for (int j = start_y; j < start_y + block_y; j++) A[i][j] += k*B[i - start_x][j - start_y];
	}
	return A;
}
//将A中分块矩阵替换为k*B
vector<vector<double>> Matrix_block_replace(vector<vector<double>> A, int start_x, int start_y, int block_x, int block_y, vector<vector<double>>B, int k) {
	for (int i = start_x; i < start_x + block_x; i++) {
		for (int j = start_y; j < start_y + block_y; j++) A[i][j] = k * B[i - start_x][j - start_y];
	}
	return A;
}


//读入m * n矩阵（）
vector<vector<double>> read_Matrix(int m,int n) {
	vector<vector<double>> A;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cin >> A[i][j];
		}
	}
	return A;
}
//求解n阶下三角方程组,前代法
vector<double> Lower_triangle(vector<vector<double>> L, vector<double> b){
	int n = L.size();
	for (int j = 0; j < n - 1; j++) {
		b[j] /= L[j][j];
		for (int k = j + 1; k <= n - 1; k++)
			b[k] -= b[j] * L[k][j];
	}
	b[n - 1] /= L[n - 1][n - 1];
	return b;
}
//求解n阶上三角方程组,回代法
vector<double> Upper_triangle(vector<vector<double>> U, vector<double> b) {
	int n = U.size();
	for (int j = n-1; j >  0; j--) {
		b[j] /= U[j][j];
		for (int k = 0; k <= j-1; k++)
			b[k] -= b[j] * U[k][j];
	}
	b[0] /= U[0][0];
	return b;
}
//输出矩阵
void cout_Matrix(vector<vector<double>> A) {
	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A[0].size(); j++) {
			printf("%6.4f",A[i][j]);
			cout << " ";
		}
			
		cout << "\n";
	}
}
//输出列向量
void cout_Lvector(vector<double> b) {
	for (int j = 0; j < b.size(); j++) {
		cout << b[j];
		cout << "\n";
	}
	cout << "\n";
}
//输出行向量
void cout_Hvector(vector<double> b) {
	for (int j = 0; j < b.size(); j++) {
		if ((j + 1) % 10 == 1) cout << "\n";
		printf("%7.4f", b[j]);
		cout << " ";
	}
}
//交换m * n阶方阵第k行和第p行
void swap_row(vector<vector<double>> &A ,int k, int p) {
	double b;
	for (int i = 0; i < A[0].size(); i++) {
		b = A[k][i];
		A[k][i] = A[p][i];
		A[p][i] = b;
	}
}
//交换m * n阶方阵第k列和第p列
void swap_col(vector<vector<double>>& A,int k, int p) {
	double b;
	for (int i = 0; i < A[0].size(); i++) {
		b = A[i][k];
		A[i][k] = A[i][p];
		A[i][p] = b;
	}
}
//找某个子矩阵（jt-->(m-1)(n-1)）里绝对值的最大元素位置
void max_subMatrix(vector<vector<double>> A, int j,int t, int m, int n, int &p, int &q) {
	for (int i = j; i <= m - 1; i++) {
		for (int k = t; k <= n - 1; k++) {
			if (abs(A[i][k]) > abs(A[p][q])) {
				p = i;
				q = k;
			}
		}
	}
}
//创建b阶单位矩阵
vector<vector<double>> I_matrix(int b) {
	vector<double> m(b);
	vector<vector<double>> n(b,m);
	for (int i = 0; i < b; i++)
		n[i][i] = 1;
	return n;
}
//交换向量b的第k,p个元素
void swap_vector(vector<double> &b, int k, int p) {
	double h;
	h = b[k]; b[k] = b[p]; b[p] = h;
}
//全主元Gauss消去法
vector<double> Gauss_all(vector<vector<double>> A,  vector<double> b) {
	int p,q,i,n = A.size();
	vector<vector<double>> L, U;
	L = make_zero_matrix(n, n);
	U = make_zero_matrix(n, n);
	vector<int> u(n), v(n);
	for (i = 0; i < n; i++) {
		u[i] = i; v[i] = i;
	}
	for (int k = 0; k < n - 1; k++) {
		p = k; q = k;
		max_subMatrix(A, k, k, n, n, p, q);
		swap_row(A, k, p);
		swap_col(A, k, q);
		u[k] = p; v[k] = q;
		if (A[k][k] != 0) {
			for (i = k + 1; i < n; i++) {
				A[i][k] /= A[k][k];
				for (int j = k + 1; j < n; j++)
					A[i][j] -= A[i][k] * A[k][j];
			}
		}
		else {
			cout << "矩阵奇异！停止运算！";
			return b;
		}
		swap_vector(b, k, p);
	}
	getLU(A,  L, U);
	b = Lower_triangle( L, b);
	b = Upper_triangle( U, b);
	//cout_Matrix(L, n, n);
	//cout_Matrix(U, n, n);
	for (int k = n-1; k >= 0; k--)
		swap_vector(b, k, v[k]);
	return b;
}
////列主元Gauss消去法
vector<double> Gauss_col(vector<vector<double>> A,  vector<double> b) {
	int p, q, i, n = A.size();
	vector<vector<double>> L, U;
	L = make_zero_matrix(n, n);
	U = make_zero_matrix(n, n);
	vector<int> u(n);
	for (i = 0; i < n; i++) {
		u[i] = i;
	}
	for (int k = 0; k < n - 1; k++) {
		p = k; q = k;
		max_subMatrix(A, k, k, n, k + 1, p, q);
		swap_row(A, k, p);
		u[k] = p;
		if (A[k][k] != 0) {
			for (i = k + 1; i < n; i++) {
				A[i][k] /= A[k][k];
				for (int j = k + 1; j < n; j++)
					A[i][j] -= A[i][k] * A[k][j];
			}
		}
		else {
			cout << "矩阵奇异！停止运算！";
			return b;
		}
		swap_vector(b, k, p);
	}
	getLU(A, L, U);

	b = Lower_triangle( L, b);
	b = Upper_triangle( U, b);
	//for (int k = n-1; k >=0; k--)
	//	swap_vector(b, k, u[k]);
	//cout_Matrix(L, n, n);
	//cout_Matrix(U, n, n);
	return b;
}

////列主元三角分解求L，U
//void Gauss_col_LU(vector<vector<double>> A, const int& n, vector<double> &L, vector<double> &U) {
//	int p, q, i;
//	vector<int> u(n);
//	for (i = 0; i < n; i++) {
//		u[i] = i;
//	}
//	for (int k = 0; k < n - 1; k++) {
//		p = k; q = k;
//		max_subMatrix(A, k, k, n, k + 1, p, q);
//		swap_row(A, n, n, k, p);
//		u[k] = p;
//		if (A[k][k] != 0) {
//			for (i = k + 1; i < n; i++) {
//				A[i][k] /= A[k][k];
//				for (int j = k + 1; j < n; j++)
//					A[i][j] -= A[i][k] * A[k][j];
//			}
//		}
//		else {
//			cout << "矩阵奇异！停止运算！";
//			return b;
//		}
//		swap_vector(b, k, p);
//	}
//	getLU(A, n, L, U);
//}
//平方根消去法
vector<double> Square_LLT(vector<vector<double>> A, vector<double> b) {
	vector<vector<double>> L, LT;
	int n = A.size();
	for (int k = 0; k < n; k++) {
		A[k][k] = sqrt(A[k][k]);
		for (int i = k + 1; i < n; i++)
			A[i][k] /= A[k][k];
		for (int j = k + 1; j < n; j++) {
			for (int s = j; s < n; s++)
				A[s][j] -= A[s][k] * A[j][k];
		}
	}
	L = make_zero_matrix(n, n);
	LT = make_zero_matrix(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			if (i >= j) L[i][j] = A[i][j];
	}
	LT = Matrix_T(L);
	//cout_Matrix(A, n, n);
	//cout << "L" << "\n";
	//cout_Matrix(L, n, n);
	//cout << "LT" << "\n";
	//cout_Matrix(LT, n, n);
	//cout << "L * LT" << "\n";
	//cout_Matrix(Matrix_Mul(L, LT, n, n, n), n, n);
	b = Lower_triangle( L, b);
	b = Upper_triangle( LT, b);
	return b;
}
//改进的平方根消去法
vector<double> Square_LDLT(vector<vector<double>> A, vector<double> b) {
	int n = A.size();
	vector<vector<double>> L, LT, D;
	vector<double> v(n);
	double w = 0;
	L = make_zero_matrix(n, n);
	LT = make_zero_matrix(n, n);
	D = make_zero_matrix(n, n);
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < j; i++) {
			v[i] = A[j][i] * A[i][i];
		}
		for (int k = 0; k < j; k++) {
			A[j][j] -= A[j][k] * v[k];
		}
		for (int s = j + 1; s < n; s++) {
			w = 0;
			for (int k = 0; k < j; k++) {
				w += A[s][k] * v[k];
			}
			A[s][j] = (A[s][j] - w) / A[j][j];
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				D[i][j] = A[i][j];
				L[i][j] = 1;
			}
			else if(i >j) L[i][j] = A[i][j];
		}
	}
	LT = Matrix_T(L);
	b = Lower_triangle( L, b);
	//cout_Matrix(L, n, n);
	b = Upper_triangle( Matrix_Mul(D,LT), b);
	return b;
}

//向量范数(k = 1,2,Infinity(-1)
double Vector_norm( vector<double> b, double k) {
	double s=abs(b[0]);
	int n = b.size();
	if (k == -1) {
		for (int i = 0; i < n; i++) 
			if(abs(b[i]) > s) s = abs(b[i]);
		return s;
	}
	else {
		for (int i = 0; i < n; i++)
			s += pow(abs(b[i]), k);
		s = pow(s, (1.0 / k));
		return s;
	}
}
//矩阵范数(k = 1,2,Infinity(-1)
double Matrix_norm(vector<vector<double>> A, int k) {
	double s=0,p=0;
	int m = A.size(), n = A[0].size();
	if (k == -1) {
		p = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++)
				p += abs(A[i][j]);
			if (p > s)s = p;
			p = 0;
		}
		return s;
	}//行和范数
	if (k == 1) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++)
				p += abs(A[i][j]);
			if (p > s)s = p;
			p = 0;
		}
		return s;
	}//列和范数
}
//找矩阵绝对值最大的元素
double Matrix_absmax(vector<vector<double>> A) {
	int m = A.size(), n = A[0].size();
	double s=A[0][0];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (abs(A[i][j]) > abs(s))s = A[i][j];
		}
	}
	return s;
}
//找矩阵最大的元素
double Matrix_max(vector<vector<double>> A) {
	int m = A.size(), n = A[0].size();
	double s = A[0][0];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (A[i][j] > s)s = A[i][j];
		}
	}
	return s;
}
//判断矩阵元素符号
vector<vector<double>> sign_Matrix(vector<vector<double>> w) {
	int m = w.size(), n = w[0].size();
	vector<vector<double>> s;
	s = make_zero_matrix(m, n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (w[i][j] > 0)s[i][j] = 1;
			else if (w[i][j] < 0)s[i][j] = -1;
		}
	}
	return s;
}
//优化法估计矩阵A逆的1范数(默认m = n,A为方阵）
double Matrix_norm_1_pro(vector<vector<double>> A) {
	int k = 1, j = 0, n = A.size();
	double f;
	vector<vector<double>> z,x,w,v,t;
	vector<double> x1(n, 0),w1(n,0),v1(n,0),z1(n,0);
	z = make_zero_matrix(n, 1);
	v = make_zero_matrix(n, 1);
	x = make_zero_matrix(n, 1);
	w = make_zero_matrix(n, 1);
	t = make_zero_matrix(1, 1);
	for (int i = 0; i < n; i++) {
		x1[i] = 1.0 / n;
		x[i][0] = 1.0 / n;
	}
	//cout_Matrix(x, n, 1);
	while (k == 1){
		j = 0;
		//B = A^-T;
		//w = Matrix_Mul(B, x, m, n, 1);
		w1 = Gauss_col(Matrix_T(A), x1);
		for (int i = 0; i < n; i++) {
			w[i][0] = w1[i];
		}
		v = sign_Matrix(w);
		//cout_Matrix(v, n, 1);
		for (int i = 0; i < n; i++) {
			v1[i] = v[i][0];
		}
		//z = Matrix_Mul(Matrix_T(B, m, n), v, n, m, 1);
		z1 = Gauss_col(A, v1);
		for (int i = 0; i < n; i++) {
			z[i][0] = z1[i];
		}
		//cout_Matrix(z, n, 1);
		t = Matrix_Mul(Matrix_T(z), x);
		if (Matrix_norm(z, -1) <= t[0][0]) {
			f = Matrix_norm(w, 1);
			k = 0;
		}
		else {
			for (int i = 0; i < n; i++) {
				if (abs(z[i][0]) == Matrix_norm(z, -1)) {
					j = i;
					i = n;
				}
			}
			x = make_zero_matrix(n, 1);
			x[j][0] = 1;
			for (int i = 0; i < n; i++)
				x1[i] = x[i][0];
			//cout_Matrix(x, n, 1);
			k = 1;
		}
	}
	return f;
}

//求A的条件数
double K(vector<vector<double>> A, int k) {
	int n = A.size();
	double t;
	if (k == -1) {
		t = Matrix_norm(A, -1) * Matrix_norm_1_pro(A);
		return t;
	}
}

//计算Householder变换  Hx = ae1
void House(vector<double> x, vector<double> &v, double &beta) {
	int n = x.size();
	double t, s = 0, a = 0;
	t = Vector_norm(x,-1);
	for (int i = 0; i < n; i++)
		x[i] /=  t;
	for (int i = 1; i < n; i++) {
		s += x[i] * x[i];
		v[i] = x[i];
	}
	//cout_Hvector(v, n);
	if (s == 0) beta = 0;
	else {
		a = sqrt(x[0] * x[0] + s);
		if (x[0] <= 0) v[0] = x[0] - a;
		else v[0] = -s / (x[0] + a);
		/*cout << v[0] << endl;*/
		beta = 2 * v[0] * v[0] / (s + v[0] * v[0]);
		for (int k = 1; k < n; k++)
			v[k] = v[k] / v[0];
		v[0] = 1;
		//cout_Hvector(v, n);
	}
}

//计算givens变换,返回cos和sin的值
void givens(double &c, double &s, double a, double b) {
	double t;
	if (b == 0) {
		c = 1;
		s = 0;
	}
	else {
		if (abs(b) > abs(a)) {
			t = a / b;
			s = 1 / sqrt(1 + t * t);
			c = s * t;
		}
		else {
			t = b / a;
			c = 1 / sqrt(1 + t * t);
			s = c * t;
		}
	}
}



//计算A的QR分解,其中B里存储的是H_k的v_k,d中存储的是beta_k
vector<vector<double>> Q_R(vector<vector<double>> A,vector<double>& d) {
	int m = A.size(), n = A[0].size();
	double beta;
	for (int j = 0; j < n; j++) {
		if (j < m) {
			vector<double> v(m - j,0),x(m - j, 0);
			vector<vector<double>> A2;
			A2 = make_zero_matrix(m - j , n - j );
			A2 = Matrix_block(A, j, j, m - j, n - j);
			for (int i = j; i < m; i++)
				x[i - j] = A[i][j];
			House(x, v, beta );
			//cout_Hvector(v, m - j);
			//cout << beta << endl;
			d[j] = beta;
			//cout_Hvector(v, m - j);
			A2 = Matrix_Subtracition(A2, Matrix_Mul_num(LVector_Mul_HVector(v, Vector_Mul_Matrix(v, A2)), beta));
			/*cout_Matrix(A2, m - j, n - j);
			cout << endl;*/
			for (int i = j; i < m; i++) {
				for (int k = j; k < n; k++)
					A[i][k] = A2[i - j][k - j];
			}
			/*cout_Matrix(A, m , n );
			cout << endl;*/
			for (int i = j + 1; i < m; i++) {
				A[i][j] = v[i - j];
			}
			/*cout_Matrix(A, m , n );
			cout << endl;*/
		}
	}
	return A;
}
//用QR分解解方程
vector<double> Q_R_solve(vector<vector<double>> A, vector<double> b) {
	int m = A.size(), n = A[0].size();
	vector<vector<double>> B,H,v,R;
	vector<double> d(m,0);
	B = make_zero_matrix(m, n);
	R = make_zero_matrix(n, n);
	H = make_zero_matrix(m, m);
	B = Q_R(A, d);
	v = make_zero_matrix(m, 1);
	for (int i = 0; i < n; i++) {
		for (int k = 0; k < m; k++) {
			if (k < i) v[k][0] = 0;
			else if (k == i) v[k][0] = 1;
			else v[k][0] = B[k][i];
		}
		//cout_Matrix(v, m, 1);
		H = Matrix_Subtracition(I_matrix(m), Matrix_Mul_num(Matrix_Mul(v, Matrix_T(v)), d[i]));
		b = Matrix_Mul_Vector(Matrix_T(H), b);
	}
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++)
			R[i][j] = B[i][j];
	}
	//cout_Matrix(B, m, n);
	//cout_Hvector(d, m);
	b = Upper_triangle( R, b);
	return b;
}

//Jacobi迭代法解方程,迭代一次
vector<double> Jacobi_Iterative_method(vector<vector<double>> A, vector<double> b, vector<double> x) {
	int n = A.size();
	vector<vector<double>> D, L, U, B;
	D = make_zero_matrix(n, n);
	L = make_zero_matrix(n, n);
	U = make_zero_matrix(n, n);
	B = make_zero_matrix(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i > j)L[i][j] = -A[i][j];
			else if (i == j)D[i][j] = A[i][j];
			else U[i][j] = -A[i][j];
		}
	}
	B =  Matrix_add(L, U);
	x = Lower_triangle(D, Vector_Add(Matrix_Mul_Vector(B, x), b));
	return x;

}

//G-S迭代法解方程,迭代一次
vector<double> G_S_Iterative_method(vector<vector<double>> A, vector<double> b, vector<double> x) {
	int n = A.size();
	vector<vector<double>> D, L, U, B, C;
	D = make_zero_matrix(n, n);
	L = make_zero_matrix(n, n);
	U = make_zero_matrix(n, n);
	B = make_zero_matrix(n, n);
	C = make_zero_matrix(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i > j)L[i][j] = -A[i][j];
			else if (i == j)D[i][j] = A[i][j];
			else U[i][j] = -A[i][j];
		}
	}
	C = Matrix_Subtracition(D, L);
	x = Lower_triangle( C, Vector_Add(Matrix_Mul_Vector(U, x), b));
	return x;
}

//SOR迭代法解方程,迭代一次
vector<double> SOR_Iterative_method(vector<vector<double>> A, vector<double> b, vector<double> x,double w) {
	int n = A.size();
	vector<vector<double>> D, L, U, B, C,E;
	D = make_zero_matrix(n, n);
	L = make_zero_matrix(n, n);
	U = make_zero_matrix(n, n);
	B = make_zero_matrix(n, n);
	C = make_zero_matrix(n, n);
	E = make_zero_matrix(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i > j)L[i][j] = -A[i][j];
			else if (i == j)D[i][j] = A[i][j];
			else U[i][j] = -A[i][j];
		}
	}
	C = Matrix_Subtracition(D, Matrix_Mul_num(L, w));
	E = Matrix_add(Matrix_Mul_num(D, 1 - w), Matrix_Mul_num(U, w));
	x = Lower_triangle( C, Vector_Add(Matrix_Mul_Vector(E, x), Vector_Mul_num(b, w)));
	return x;
}

//Jacobi迭代法解方程，二维偏微分方程
vector<vector<double>> Jacobi_Iterative_method2(vector<vector<double>> U, vector<vector<double>> F, vector<vector<double>> G) {
	vector<vector<double>> U0;
	int N = U.size() - 1;
	double h = 1 / N;
	U0 = make_zero_matrix(N+1,N+1);
	U0 = U;
	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++)
			U0[i][j] = (F[i][j] + U[i - 1][j] + U[i][j - 1] + U[i + 1][j] + U[i][j + 1]) / (4 + G[i][j]);
	}
	return U0;
}

//G-S迭代法解方程,二维偏微分方程
vector<vector<double>> G_S_Iterative_method2(vector<vector<double>> U, vector<vector<double>> F, vector<vector<double>> G) {
	int N = U.size() - 1;
	double h = 1 / N;
	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++)
			U[i][j] = (F[i][j] + U[i - 1][j] + U[i][j - 1] + U[i + 1][j] + U[i][j + 1]) / (4 + G[i][j]);
	}
	return U;
}

//SOR迭代法解方程,二维偏微分方程
vector<vector<double>> SOR_Iterative_method2(vector<vector<double>> U, vector<vector<double>> F, vector<vector<double>> G, double w) {
	int N = U.size() - 1;
	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++)
			U[i][j] = (w * (F[i][j] + U[i - 1][j] + U[i][j - 1] + U[i + 1][j] + U[i][j + 1]) - (w - 1) * (4 +  G[i][j])*U[i][j]) / (4 + G[i][j]);
	}
	return U;
}

//共轭梯度法,迭代一次
void conjugate_grad_method_1(vector<vector<double>> A, vector<double> b, vector<double> &x, vector<double> &p, vector<double> &r) {
	int n = A.size();
	double i, j;
	i = HVector_Mul_LVector(r, p) / HVector_Mul_LVector(p, Matrix_Mul_Vector(A, p));
	x = Vector_Add(x, Vector_Mul_num(p, i));
	r = Vector_Subtracition(b, Matrix_Mul_Vector(A, x));
	j = 0 - HVector_Mul_LVector(r, Matrix_Mul_Vector(A, p)) / HVector_Mul_LVector(p, Matrix_Mul_Vector(A, p));
	p = Vector_Add(r, Vector_Mul_num(p, j));

}

//共轭梯度法求解,迭代停止条件为||x1-x2||<w,返回迭代次数k
vector<double> conjugate_grad_method(vector<vector<double>> A, vector<double> b, double w,int &k) {
	int n = A.size();
	vector<double>x(n, 0.5), x1(n), x2(n), r(n), p(n);
	double  a, beta;
	r = Vector_Subtracition(b, Matrix_Mul_Vector(A, x));                            //r0 = b-Ax0
	p = r;                                                                          //p0 = r0
	a = HVector_Mul_LVector(r, r) / HVector_Mul_LVector(p, Matrix_Mul_Vector(A, p));//a0 = r0Tr0/p0TAp0
	x1 = Vector_Add(x, Vector_Mul_num(p, a));                                        //x1 = x0 + a0p0
	r = Vector_Subtracition(b, Matrix_Mul_Vector(A, x1));                           //r1 = b-Ax1
	beta = -HVector_Mul_LVector(r, Matrix_Mul_Vector(A, p)) / HVector_Mul_LVector(p, Matrix_Mul_Vector(A, p));//beta0=r1TAp0/p0TAp0
	p = Vector_Add(r, Vector_Mul_num(p, beta));                                        //p1 = r1 +beta0*p0
	x2 = x;
	while ((Vector_norm(Vector_Subtracition(x1, x2), -1)) > pow(10, -7)) {
		x2 = x1;
		conjugate_grad_method_1(A, b, x1, p, r);
		k++;
	}
	return x1;
}
//幂法迭代一次
double Power_method_1(vector<vector<double>> A, vector<double>& u) {
	double lanta;
	vector<double> y(u.size());
	y = Matrix_Mul_Vector(A, u);
	lanta = y[0];
	for (int i = 1; i < u.size(); i++) 
		if (abs(y[i]) > abs(lanta)) lanta = y[i];
	u = Vector_Mul_num(y, 1.0 / lanta);
	return lanta;
}
//幂法求模最大特征值和特征向量
double Power_method(vector<vector<double>> A, vector<double>& u, double w,int &k) {
	double lanta;
	k = 0;
	vector<double> y(u.size()),u0(u.size());
	u0 = u;
	lanta = Power_method_1(A, u);
	while (Vector_norm(Vector_Subtracition(u0, u), -1) > w) {
		k++;
		u0 = u;
		lanta = Power_method_1(A, u);
		//cout << lanta << endl;
	}
	return lanta;
}
//Household变换法求上Hessenberg分解（算法6.4.1）次对角线要非零
vector<vector<double>> up_Hessenberg_with_household(vector<vector<double>> A) {
	int n = A.size();
	double beta = 0;
	for (int k = 0; k <= n - 3; k++) {
		vector<double> v(n - k - 1), x(n - k - 1);
		for (int j = k + 1; j < n; j++) x[j - k - 1] = A[j][k];
		House(x, v, beta);
		vector<vector<double>> I = I_matrix(n - k - 1);
		A = Matrix_block_mul(A, k + 1, k, n - k - 1, n - k, Matrix_Subtracition(I, Matrix_Mul_num(LVector_Mul_HVector(v, v), beta)), -1);
		A = Matrix_block_mul(A, 0, k + 1, n, n - k - 1, Matrix_Subtracition(I, Matrix_Mul_num(LVector_Mul_HVector(v, v), beta)), 1);
	}
	return A;
}
//算法6.4.2双重步位移的QR迭代,返回变换后的H,保存P
vector<vector<double>> QR_iteration(vector<vector<double>>H, vector<vector<double>> &P) {
	int n = H.size(), m = n - 1, q, r;
	double s, t, x, y, z, beta = 0;
	vector<vector<double>> I3 = I_matrix(3),P1,H1=I3,H0=I_matrix(2);
	s = H[m - 1][m - 1] + H[n - 1][n - 1];
	t = H[m - 1][m - 1] * H[n - 1][n - 1] - H[m - 1][n - 1] * H[n - 1][m - 1];
	x = H[0][0] * H[0][0] + H[0][1] * H[1][0] - s * H[0][0] + t;
	y = H[1][0] * (H[0][0] + H[1][1] - s);
	z = H[1][0] * H[2][1];
	for (int k = -1; k < n - 3; k++) {
		vector<double> w = { x,y,z }, v(3);
		House(w, v, beta);
		P1 = I_matrix(n);
		H1 = Matrix_Subtracition(I3, Matrix_Mul_num(LVector_Mul_HVector(v, v), beta));
		for (int i = k + 1; i < k + 4; i++) {
			for (int j = k + 1; j < k + 4; j++) {
				P1[i][j] = H1[i - k - 1][j - k - 1];
			}
		}
		P = Matrix_Mul(P, P1);
		q = 0; if (k > 0)q = k;
		H = Matrix_block_mul(H, k + 1, q, 3, n - q, H1, -1);
		r = k + 4; if (r > n - 1) r = n - 1;
		H = Matrix_block_mul(H, 0, k+1, r+1, 3, H1, 1);
		x = H[k + 2][k + 1];
		y = H[k + 3][k + 1];
		if (k < n - 4) z = H[k + 4][k + 1];
	}
	vector<double> w0 = { x,y }, v0(2);
	vector<vector<double>> I2 = I_matrix(2);
	House(w0, v0, beta);
	P1 = I_matrix(n);
	H0 = Matrix_Subtracition(I2, Matrix_Mul_num(LVector_Mul_HVector(v0, v0), beta));
	for (int i = n-2; i < n; i++) {
		for (int j = n-2; j < n; j++) {
			P1[i][j] = H0[i - n + 2][j - n + 2];
		}
	}
	P = Matrix_Mul(P, P1);
	H = Matrix_block_mul(H, n-2, n-3, 2, 3, Matrix_Subtracition(I2, Matrix_Mul_num(LVector_Mul_HVector(v0, v0), beta)), -1);
	H = Matrix_block_mul(H, 0, n - 2, n, 2, Matrix_Subtracition(I2, Matrix_Mul_num(LVector_Mul_HVector(v0, v0), beta)), 1);
	return H;
}
//算法6.4.3确定最大的m和l
void m_l(vector<vector<double>> H, int &m, int &l) {
	int n = H.size(),k;
	//m
	for (m = 1; m < n; m++) {
		if (H[n - m][n - m - 1] == 0)continue;
		else {
			if (m == n - 1) {
				if (pow(H[n - m - 1][n - m - 1] + H[n - m][n - m], 2) - 4 * (H[n - m - 1][n - m - 1] * H[n - m][n - m] - H[n - m - 1][n - m] * H[n - m][n - m - 1]) < 0) {
					m = n;
					break;
				}
			}
			else {
				if (H[n - m - 1][n - m - 2] != 0) {
					m = m - 1;
					break;
				}
				else {
					if (pow(H[n - m - 1][n - m - 1] + H[n - m][n - m], 2) - 4 * (H[n - m - 1][n - m - 1] * H[n - m][n - m] - H[n - m - 1][n - m] * H[n - m][n - m - 1]) < 0)continue;
					else {
						m = m - 1;
						break;
					}
				}
			}
		}
	}
	//l
	l = 0;
	while ((H[l+1][l]==0) && (l<n-m+1))
	{
		l++;
	}
	
	
	
}


//隐式QR算法求解A的特征值,机器精度为u
vector<vector<double>> Implicit_QR_algo(vector<vector<double>> A, double u,int &k) {
	int n = A.size(),l=0,m=0;
	k = 0;
	vector<vector<double>> H,H22,H12,H23,P;
	double t, sin, cos, a, b, c, d;
	H = make_zero_matrix(n, n);
	H = up_Hessenberg_with_household(A);
	for (int i = 1; i < n; i++) if (fabs(H[i][i - 1]) <= (fabs(H[i][i]) + fabs(H[i - 1][i - 1])) * u) H[i][i - 1] = 0;
	m_l(H, m, l);
	//cout <<"986" << m <<" " << n << endl;
	while (m != n)
	{
		//cout << "989" << m << " " << n << endl;
		//cout << "990" << m <<" " << n << endl;
		if (n - l - m > 2) {
			H22 = make_zero_matrix(n - l - m, n - l - m);
			P = I_matrix(n - l - m);
			H22 = Matrix_block(H, l, l, n - l - m, n - l - m);
			H22 = QR_iteration(H22, P);
		}
		else {
			P = make_zero_matrix(2, 2);
			a = H22[0][0]; b = H22[0][1]; c = H22[1][0]; d = H22[1][1];
			t = ((a - d) + sqrt(pow(a - d, 2) + 4 * b * c)) / (2 * b);
			sin = t / sqrt(1 + t * t);
			cos = 1 / sqrt(1 + t * t);
			P[0][0] = P[1][1] = cos;
			P[0][1] = sin;
			P[1][0] = -sin;
		}
		if (l > 0) {
			H = Matrix_block_mul(H, 0, l, l, n - l - m, P, 1);//H12=H12*P
		}
		H = Matrix_block_mul(H, l, l, n - l - m, n - l - m, Matrix_T(P), -1);//H22=PT*H22
		H = Matrix_block_mul(H, l, l, n - l - m, n - l - m, P, 1);//H22=H22*P
		H = Matrix_block_mul(H, l, n - m, n - l - m, m, Matrix_T(P), -1);//H23=PT*H23
		for (int i = 1; i < n; i++) if (fabs(H[i][i - 1]) <= (fabs(H[i][i]) + fabs(H[i - 1][i - 1])) * u) H[i][i - 1] = 0;
		m_l(H, m, l);
		k++;
		//cout_Matrix(H);
		//cout << m << l << A[1][0] << endl;
		//cout << endl;
	}
	return H;
}
//输出拟对角阵的特征值
void cout_diag_characteristic_value(vector<vector<double>> A) {
	int n = A.size();
	double a, b, c, d, p, q, s, w=0;
	for (int i = 0; i < n - 1; i++) {
		if (A[i + 1][i] < pow(10,-5)) cout << A[i][i] << endl;
		else {
			a = A[i][i]; b = A[i][i + 1]; c = A[i + 1][i]; d = A[i + 1][i + 1];
			s = pow(a + d, 2) - 4 * (a * d - b * c);
			p = (a + d) / 2.0;
			q = sqrt(-s) / 2.0;
			cout << p << "+" << q << "i" << "  " << "\n" << p << "-" << q << "i" << endl;
			if (i == n - 2) w = 1;
			else i++;
		}
	}
	if (w != 1)cout << A[n - 1][n - 1];
}

//乘上givens变换,k=1左乘,k=-1右乘
vector<vector<double>> givens_mul(vector<vector<double>> A, double cos, double sin, int p, int q, int k) {
	int n = A.size();
	vector<vector<double>> B;
	B = make_zero_matrix(n,n);
	B = A;
	if (k == 1) {
		for (int i = 0; i < n; i++) {
			B[p][i] = cos * A[p][i] + sin * A[q][i];
			B[q][i] = -sin * A[p][i] + cos * A[q][i];
		}
	}
	else {
		for (int i = 0; i < n; i++) {
			B[i][p] = cos * A[i][p] - sin * A[i][q];
			B[i][q] = sin * A[i][p] + cos * A[i][q];
		}
	}
	return B;
}
//计算A的非对角范数
double Non_diagonal_norm(vector<vector<double>> A) {
	int n = A.size();
	double E = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) E = E;
			else E = E + A[i][j] * A[i][j];
		}
	}
	return E;
}

//过关Jacobi法求实对称矩阵特征值和特征向量，保存在A和J中
vector<vector<double>> Pass_Jacobi(vector<vector<double>>A, vector<vector<double>>& J, int& k) {
	int n = A.size(), p = 0, q = 0, a = 0;
	double cos, sin, t, m, E = Non_diagonal_norm(A), h = E;
	k = 0;
	while (E > pow(10, -12)) {
		a = 0;
		//cout << "E=" << E << endl;
		for (int q = 1; q < n; q++) {
			for (int p = 0; p < q; p++) {
				if (fabs(A[q][p]) > h) {
						a++;
						m = (A[q][q] - A[p][p]) / (2 * A[p][q]);
						t = sgn(m) /  (fabs(m) + sqrt(1 + m * m));
						cos = 1 / sqrt(1 + t * t);
						sin = t * cos;//确定cos和sin
						J = givens_mul(J, cos, sin, p, q, -1);//J=J1J2J3J4....
						A = givens_mul(givens_mul(A, cos, -sin, p, q, 1), cos, sin, p, q, -1);//A=JT*A*J
						k++;
						//cout_Matrix(A);
						//cout << endl;
				}
			}
		}//确定p,q
		if (a == 0) {
			h /= n;
			E = Non_diagonal_norm(A);
		}
	}
	return A;
}
//符号判断函数
int sgn(double m) {
	if (m >= 0)return 1;
	else return -1;
}
 //符号判断函数,  >0为1   <0为-1 ,0为0
int sign(double m) {
	if (m > 0)return 1;
	else if (m < 0)return -1;
	else return 0;
}

//计算变号数sn(u)
int variants_num(vector<vector<double>>A, double u) {
	int n = A.size(), s = 0;
	vector<double>x(n), y(n);
	x[0] = A[0][0];
	y[0] = 0;
	for (int i = 1; i < n; i++) {
		x[i] = A[i][i];
		y[i] = A[i][i - 1];
	}//初始化x,y
	double q = x[0] - u;
	for (int k = 0; k < n; k++) {
		if (q < 0) s++;
		if (k < n - 1) {
			if (q == 0)q = fabs(y[k + 1]) * pow(10, -10);
			q = x[k + 1] - u - y[k + 1] * y[k + 1] / q;
		}
	}
	return s;
}
//用二分法求一个实对称三对角阵A中第m大的特征值
double dichotomy(vector<vector<double>>A, int m,int &k) {
	int n = A.size();
	k = 0;
	double l = -Matrix_norm(A, -1), u = Matrix_norm(A, -1),r;
	while (u - l > pow(10, -10)) {
		r = (u + l) / 2;
		if (variants_num(A, r) >= m)u = r;
		else l = r;
		k++;
	}
	return (l + u) / 2;
}
//反幂法求lenta的特征向量
vector<double> Inverse_Power_method(vector<vector<double>>A, double lanta) {
	int n = A.size();
	vector<vector<double>>B, I;
	B = I = I_matrix(n);
	B = Matrix_Subtracition(A, Matrix_Mul_num(I, lanta));
	vector<double> z0(n, 1), z1(n, 0),v(n);
	while (Vector_norm(Vector_Subtracition(z0, z1), -1)>pow(10,-6)) {
		z1 = z0;
		v = Gauss_col(B, z0);
		double s = v[0];
		for (int i = 1; i < v.size(); i++)
			if (abs(v[i]) > abs(s)) s = v[i];
		z0 = Vector_Mul_num(v, 1 / s);
		//cout_Hvector(z0);
	}
	return z0;
}
//二对角化
vector<vector<double>> two_diagonal(vector<vector<double>>A, vector<vector<double>>& U, vector<vector<double>>& V) {
	int m = A.size(), n = A[0].size();
	double beta = 0;
	vector<double> b(n), c(n-1);
	vector<vector<double>> U1,V1,B,I;
	V = I_matrix(n);
	U = I_matrix(m);
	for (int k = 0; k < n; k++) {
		vector<double> v(m - k), u(n - k), x(m - k);
		for (int i = k; i < m; i++) x[i - k] = A[i][k];
		House(x, v, beta);
		u = Vector_Mul_Matrix(Vector_Mul_num(v, beta), Matrix_block(A, k, k, m - k, n - k));
		A = Matrix_block_add(A, k, k, m - k, n - k, LVector_Mul_HVector(v, u), -1);
		for (int i = k + 1; i < m; i++) A[i][k] = v[i - k];
		b[k] = beta;
		u = Vector_Mul_Matrix(Vector_Mul_num(v, beta), Matrix_block(U, k, 0, m - k, m));
		U = Matrix_block_add(U, k, 0, m - k, m, LVector_Mul_HVector(v, u), -1);
		//cout_Matrix(U); cout << "\n" << endl;
		if (k < n - 2) {
			vector<double> v(n - 1 - k), u(m - k), x(n - 1 - k);
			for (int i = k + 1; i < n; i++) x[i - k - 1] = A[k][i];
			House(x, v, beta);
			u = Matrix_Mul_Vector(Matrix_block(A, k, k+1, m - k, n - k - 1), Vector_Mul_num(v, beta));
			A = Matrix_block_add(A, k, k+1, m - k, n - k - 1, LVector_Mul_HVector(u, v), -1);
			for (int i = k + 2; i < n; i++) A[k][i] = v[i - k - 1];
			c[k] = beta;
			u = Matrix_Mul_Vector(Matrix_block(V, 0, k + 1, n, n - k - 1), Vector_Mul_num(v, beta));
			V = Matrix_block_add(V, 0, k + 1, n, n - k - 1, LVector_Mul_HVector(u, v), -1);
			//cout_Matrix(V); cout << "\n" << endl;
		}
	}
	c[n - 2] = beta;
	B = I_matrix(n);
	B[0][0] = A[0][0];
	for (int i = 1; i < n; i++) {
		B[i][i] = A[i][i];
		B[i - 1][i] = A[i - 1][i];
	}
	//cout_Matrix(A);
	return B;
}
//SVD迭代中求pq
void p_q(vector<vector<double>>B, int& p, int& q) {
	int n = B.size(),l;
	p = q = 0;
	for (q = 0; q < n - 1; q++) {
		//cout << "1" << endl;
		if (B[n - 2 - q][n - 1 - q] != 0)break;
	}
	if (q == n - 1)q = n;
	else {
		for (l = 1; l < n - q - 1; l++) {
			//cout << "2" << endl;
			if (B[n - q - 2 - l][n - q - 1 - l] == 0)break;
		}
		p = n - q - l - 1;
	}
}
//算法7.6.2带Wilkinson位移的SVD迭代   PT*B*Q
vector<vector<double>> Wilkinson_SVD(vector<vector<double>>B, vector<vector<double>>& P, vector<vector<double>>& Q) {
	int n = B.size();
	double alpha, deta, beta, mu, gama, y, z, c, s;
	alpha = pow(B[n - 1][n - 1], 2) + pow(B[n - 2][n - 1], 2);
	deta = (pow(B[n - 2][n - 2], 2) + pow(B[n - 3][n - 2], 2) - alpha) / 2;
	beta = B[n - 2][n - 2] * B[n - 2][n - 1];
	mu = alpha - beta * beta / (deta + sign(deta) * sqrt(deta * deta + beta * beta));
	y = pow(B[0][0], 2) - mu;
	z = B[0][0] * B[0][1];
	for (int k = 0; k < n - 1; k++) {
		givens(c, s, y, z);
		B = givens_mul(B, c, -s, k, k + 1, -1);
		Q = givens_mul(Q, c, -s, k, k + 1, -1);
		y = B[k][k];
		z = B[k + 1][k];
		givens(c, s, y, z);
		P = givens_mul(P, c, s, k, k + 1, 1);
		if (k < n - 2) {
			B = givens_mul(B, c, s, k, k + 1, 1);
			y = B[k][k + 1];
			z = B[k][k + 2];
		}
		else {
			B = givens_mul(B, c, s, k, k + 1, 1);

			//y = c * B[n - 2][n - 1] + s * B[n - 1][n - 1];
			//z = -s * B[n - 2][n - 1] + c * B[n - 1][n - 1];
			//B[n - 2][n - 1] = y;
			//B[n - 1][n - 1] = z;
		}
	}
	return B;
}
//输出B的奇异值,k>0升序，k=0直接输出，k<0降序
void cout_Singular_values(vector<vector<double>>B, int k) {
	int n = B.size();
	double t;
	if (k == 0) {
		for (int i = 0; i < n; i++) if (B[i][i] > 0)cout << B[i][i] << " ";
	}
	else if (k != 0) {
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				if (k * B[j][j] < k * B[i][i]) {
					t = B[i][i];
					B[i][i] = B[j][j];
					B[j][j] = t;
				}
			}
		}
		cout_Singular_values(B, 0);
	}
}














