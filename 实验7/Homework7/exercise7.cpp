#include <iostream>
#include <vector>
#include"Function.h"
#include"exercise7.h"
#include<time.h>
//#include<random>
#include<cstdlib>

using namespace std;

void exercise_7_1() {
	//int n = 4;
	//vector<vector<double>>A, J;
	//A = J = I_matrix(n);
	//A[0][0] = 4;
	//for (int i = 1; i < n; i++) {
	//	A[i][i - 1] = A[i - 1][i] = 1;
	//	A[i][i] = 4;
	//}//初始化A
	//cout_Matrix(A);
	//cout_Matrix(givens_mul(A, 0.6, 0.8, 2, 3, 1));

	int k;
	clock_t start, finish;
	double t;
	srand(time(0));
	for (int n = 50; n < 101; n += 10) {
		vector<vector<double>>A, J;
		A = J = I_matrix(n);
		A[0][0] = 4;
		for (int i = 1; i < n; i++) {
			A[i][i - 1] = A[i - 1][i] = 1;
			A[i][i] = 4;
		}//初始化A
		start = clock();
		A = Pass_Jacobi(A, J, k);
		finish = clock();
		t = finish - start;
		cout << "n=" << n << "  迭代次数:" << k << "  用时" << t << "ms" << endl;
		cout << "Ak=" << endl;
		cout_Matrix(A);
		cout << "\n" << "Qk=" << endl;
		cout_Matrix(J);
		if (n > 50) {
			double t;
			for (int i = 0; i < n; i++) {
				for (int j = i + 1; j < n; j++) {
					if (A[j][j] < A[i][i]) {
						t = A[i][i];
						A[i][i] = A[j][j];
						A[j][j] = t;
					}
				}
			}
			cout << endl;
			cout_diag_characteristic_value(A);
		}
			
	}
}




void exercise_7_2() {
	int k, n = 100;
	clock_t start, finish;
	double t,lanta;
	srand(time(0));
	vector<vector<double>>A, J;
	A = J = I_matrix(n);
	A[0][0] = 2;
	for (int i = 1; i < n; i++) {
		A[i][i - 1] = A[i - 1][i] = -1;
		A[i][i] = 2;
	}//初始化A
	//cout_Matrix(A);
	for (int m = 1; m <= n; m += n - 1) {
		start = clock();
		lanta = dichotomy(A, m, k);
		finish = clock();
		t = finish - start;
		if (m == 1) {
			cout << "最小特征值:" << lanta << "  迭代次数:" << k << "  用时:" << t << "ms" << endl;
			cout << "特征向量:" << endl;
			cout_Lvector(Inverse_Power_method(A, lanta));
		}
		else {
			cout << "最大特征值:" << lanta << "  迭代次数:" << k << "  用时:" << t << "ms" << endl;
			cout << "特征向量:" << endl;
			cout_Lvector(Inverse_Power_method(A, lanta));
		}
	}
}























