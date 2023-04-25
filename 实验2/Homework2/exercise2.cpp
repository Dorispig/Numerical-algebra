#include <iostream>
#include <vector>
#include"Function.h"
#include"exercise2.0.h"
#include<time.h>
using namespace std;
//ϰ��2
void exercise_2_1(void) {
	for (int m = 1; m < 21; m++) {
	double k=0;
	vector<vector<double>> A,L,U;
	A = make_zero_matrix(m, m);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++)
			A[i][j] = 1.0 /(i+j+1);
	}//Hilbert����
	//L = make_zero_matrix(m, m);
	//U = make_zero_matrix(m, m);
	//Gauss_col_LU(A, m, L, U);
	k = Matrix_norm(A, m, m, -1) * Matrix_norm_1_pro(A, m);
	cout << "��Ϊ" << m << "��Hilbert����������������Ϊ��" << k << endl;
	}
}
 
void exercise_2_2(void) {
	srand(time(0));
	double t;
	for (int n = 5; n < 31; n++) {
		t = 0;
		vector<vector<double>> A;
		vector<double> x(n, 0),b(n,0),x1(n,0);
		A = make_zero_matrix(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j || j == n - 1) A[i][j] = 1;
				else if (i > j) A[i][j] = -1;
			}
			x[i] = double(rand()) / (RAND_MAX) * 100.0;
		}//����A����,x����
		b = Matrix_Mul_Vector(A, x, n, n);
		x1 = Gauss_col(A, n, b);
		t = K(A, n, -1) * Vector_norm(Vector_Subtracition(b, Matrix_Mul_Vector(A, x1, n, n), n), n, -1) / Vector_norm(b, n, -1);
		cout << "��Ϊ" << n << "�ľ���Ϊ��" << t << "       ";
		t = Vector_norm(Vector_Subtracition(x, x1, n), n, -1) / Vector_norm(x, n, -1);
		cout << "������Ϊ��" << t<<endl;
	}
}







































