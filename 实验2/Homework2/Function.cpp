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
vector<vector<double>> Matrix_Mul(vector<vector<double>> a, vector<vector<double>> b,const int& m,const int& n,const int& s) {
	vector<vector<double>> p;
	p = make_zero_matrix(m,n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < s; j++) {
			for (int k = 0; k < n; k++) 
				p[i][j] += a[i][k] * b[k][j];
		}
	}
	return p;
}

//求矩阵转置
vector<vector<double>> Matrix_T(vector<vector<double>> a,const int& m, const int& n) {
	vector<vector<double>> b;
	b = make_zero_matrix(n,m);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			b[j][i] = a[i][j];
	}
	return b;
}

//矩阵数乘
vector<vector<double>> Matrix_Mul_num(vector<vector<double>> a, double b) {
	for (int i = 0; i < a.size(); i++) {
		for (int j = 0; j < a[0].size(); j++)
			a[i][j] *= b;
	}
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
//矩阵相减
vector<vector<double>> Matrix_Subtracition(vector<vector<double>> a, vector<vector<double>> b) {
	for (int i = 0; i < a.size(); i++) {
		for (int j = 0; j < a[0].size(); j++)
			a[i][j] -= b[i][j];
	}
	return a;
}
//向量相减
vector<double> Vector_Subtracition(vector<double> a, vector<double> b,int m) {
	vector<double> c(m);
	for (int i = 0; i < m; i++) 
			c[i] = a[i] - b[i];
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
	double n=1;
	if (A.size() == 1)
		return A[0][0];
	else {
			for (int i = 0; i < A.size() - 1; i++) {
				for (int j = i + 1; j < A.size(); j++) {
					for (int k = i; k < A.size(); k++)
						A[j][k] -= A[i][k] * A[j][i] / A[i][i];
				}
			}
			for (int i = 0; i < A.size(); i++)
				n *= A[i][i];
			return n;
	}	
}
//不选主元的Gauss消去法求L,U,直接用for循环
vector<double> Gauss_no(vector<vector<double>> A,const int &n,vector<double> b){
	vector<vector<double>>L;
	vector<vector<double>>U;
	L = make_zero_matrix(n, n);
	U = make_zero_matrix(n, n);
	for (int k = 0; k < n-1; k++) {
		for (int i = k + 1; i < n; i++) {
			A[i][k] /= A[k][k];
			for (int j = k + 1; j < n; j++)
				A[i][j] -= A[i][k] * A[k][j];
		}
	}
	getLU(A, n, L, U);
	b = Lower_triangle(n, L, b);
	b = Upper_triangle(n, U, b);
	return b;
}
//取L,U矩阵
void getLU(vector<vector<double>> A, const int& n, vector<vector<double>>& L, vector<vector<double>>& U) {
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
	for (int i = 0; i < block_x; i++) {
		for (int j = 0; j < block_y; j++) {
			B[i][j] = A[start_x + i][start_y + j];
		}
	}
	return B;
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
vector<double> Lower_triangle(int n,vector<vector<double>> L, vector<double> b){
	for (int j = 0; j < n - 1; j++) {
		b[j] /= L[j][j];
		for (int k = j + 1; k <= n - 1; k++)
			b[k] -= b[j] * L[k][j];
	}
	b[n - 1] /= L[n - 1][n - 1];
	return b;
}
//求解n阶上三角方程组,回代法
vector<double> Upper_triangle(int n,vector<vector<double>> U, vector<double> b) {
	for (int j = n-1; j >  0; j--) {
		b[j] /= U[j][j];
		for (int k = 0; k <= j-1; k++)
			b[k] -= b[j] * U[k][j];
	}
	b[0] /= U[0][0];
	return b;
}
//输出矩阵
void cout_Matrix(vector<vector<double>> A, int m, int n) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			printf("%7.2f",A[i][j]);
			cout << " ";
		}
			
		cout << "\n";
	}
}

//输出列向量
void cout_Lvector(vector<double> b, int n) {
	for (int j = 0; j < n; j++) {
		cout << b[j];
		cout << "\n";
	}
	cout << "\n";
}
//输出行向量
void cout_Hvector(vector<double> b, int n) {
	for (int j = 0; j < n; j++) {
		printf("%7.2f", b[j]);
		cout << " ";
		if ((j + 1) % 10 == 0) cout << "\n";
	}
}
//交换m * n阶方阵第k行和第p行
void swap_row(vector<vector<double>> &A, int m, int n, int k, int p) {
	double b;
	for (int i = 0; i < n; i++) {
		b = A[k][i];
		A[k][i] = A[p][i];
		A[p][i] = b;
	}
}
//交换m * n阶方阵第k列和第p列
void swap_col(vector<vector<double>>& A, int m, int n, int k, int p) {
	double b(n);
	for (int i = 0; i < n; i++) {
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
vector<vector<double>> make_I_matrix(int b) {
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
vector<double> Gauss_all(vector<vector<double>> A, const int& n, vector<double> b) {
	int p,q,i;
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
		swap_row(A, n, n, k, p);
		swap_col(A, n, n, k, q);
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
	getLU(A, n, L, U);
	b = Lower_triangle(n, L, b);
	b = Upper_triangle(n, U, b);
	//cout_Matrix(L, n, n);
	//cout_Matrix(U, n, n);
	for (int k = n-1; k >= 0; k--)
		swap_vector(b, k, v[k]);
	return b;
}
////列主元Gauss消去法
vector<double> Gauss_col(vector<vector<double>> A, const int& n, vector<double> b) {
	int p, q, i;
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
		swap_row(A, n, n, k, p);
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
	getLU(A, n, L, U);

	b = Lower_triangle(n, L, b);
	b = Upper_triangle(n, U, b);
	//for (int k = n-1; k >=0; k--)
	//	swap_vector(b, k, u[k]);
	//cout_Matrix(L, n, n);
	//cout_Matrix(U, n, n);
	return b;
}

//列主元三角分解求L，U
void Gauss_col_LU(vector<vector<double>> A, const int& n, vector<double> &L, vector<double> &U) {
	int p, q, i;
	vector<int> u(n);
	for (i = 0; i < n; i++) {
		u[i] = i;
	}
	for (int k = 0; k < n - 1; k++) {
		p = k; q = k;
		max_subMatrix(A, k, k, n, k + 1, p, q);
		swap_row(A, n, n, k, p);
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
	getLU(A, n, L, U);
}
//平方根消去法
vector<double> Square_LLT(vector<vector<double>> A, const int& n, vector<double> b) {
	vector<vector<double>> L, LT;
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
	LT = Matrix_T(L, n);
	//cout_Matrix(A, n, n);
	//cout << "L" << "\n";
	//cout_Matrix(L, n, n);
	//cout << "LT" << "\n";
	//cout_Matrix(LT, n, n);
	//cout << "L * LT" << "\n";
	//cout_Matrix(Matrix_Mul(L, LT, n, n, n), n, n);
	b = Lower_triangle(n, L, b);
	b = Upper_triangle(n, LT, b);
	return b;
}
//改进的平方根消去法
vector<double> Square_LDLT(vector<vector<double>> A, const int& n, vector<double> b) {
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
	LT = Matrix_T(L, n);
	b = Lower_triangle(n, L, b);
	//cout_Matrix(L, n, n);
	b = Upper_triangle(n, Matrix_Mul(D,LT,n,n,n), b);
	return b;
}

//向量范数
double Vector_norm( vector<double> b, int n, double k) {
	double s=0;
	for (int i = 0; i < n; i++)
		s += pow(b[i], k);
	s = pow(s, (1.0 / k));
	return s;
}
//矩阵范数(k = 1,2,Infinity(-1)
double Matrix_norm(vector<vector<double>> A, int m, int n, int k) {
	double s=0,p=0;
	if (k == -1) {
		p = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++)
				p += abs(A[i][j]);
			if (p > s)s = p;
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
//判断矩阵元素符号
vector<vector<double>> sign_Matrix(vector<vector<double>> w, int m, int n) {
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
//优化法估计矩阵A逆的1范数(默认m = n）
double Matrix_norm_1_pro(vector<vector<double>> A,int m,int n) {
	int k = 1,j=0;
	double f;
	vector<vector<double>>  w, v, z;
	vector<double> x(n, 0), v1(m, 0);
	w = make_zero_matrix(m, 1);
	v = make_zero_matrix(m, 1);
	z = make_zero_matrix(n, 1);
	for (int i = 0; i < n; i++) 
		x[i] = 1.0 / (i + 1);
	do while (k == 1) {
		j = 0;
		//B = A^-T;
		//w = Matrix_Mul(B, x, m, n, 1);
		w = Gauss_col(Matrix_T(A, m, n), n, x);
		v = sign_Matrix(w, m, 1);
		for (int i = 0; i < m; i++) {
			v1[i] = v[i][0];
		}
		//z = Matrix_Mul(Matrix_T(B, m, n), v, n, m, 1);
		z = Gauss_col(A, n, v);
		if (Matrix_norm(z, n, 1, -1) <= Matrix_Mul(Matrix_T(z, n, 1), x, 1, n, 1)[0][0]) {
			f = Matrix_norm(z, n, 1, 1);
			k = 0;
		}
		else {
			for (int i = 0; i < n; i++) {
				if (z[i][0] == Matrix_norm(z, n, 1, -1)) {
					j = i;
					i = n;
				}
			}
			x = make_zero_matrix(n, 1);
			x[j][0] = 1;
			k = 1;
		}
	}
	return f;
}













