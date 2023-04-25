#include <iostream>
#include <vector>
#include"Function.h"
#include"exercise6.h"
#include<time.h>
//#include<random>
#include<cstdlib>

using namespace std;


void exercise_6_1_1() {
	cout << "x^3+x^2-5x+3=0的模最大根为：  ";
	vector<vector<double>> A = { {-1,5,-3},{1,0,0},{0,1,0} };
	vector<double> u(3, 0);
	int k;
	clock_t start, finish;
	double t;
	srand(time(0));
	u[0] = 1;
	double lanta;
	start = clock();
	lanta = Power_method(A, u, pow(10, -6),k);
	finish = clock();
	t = finish - start;
	cout << lanta << endl;
	cout << "迭代次数为：" << k << "   运行时间为：" << t << "ms" << endl << endl;
}
void exercise_6_1_2() {
	cout << "x^3-3x-1=0的模最大根为：  ";
	vector<vector<double>> A = { {0,3,1},{1,0,0},{0,1,0} };
	vector<double> u(3, 1);
	int k;
	clock_t start, finish;
	double t;
	srand(time(0));
	double lanta;
	start = clock();
	lanta = Power_method(A, u, pow(10, -4),k);
	finish = clock();
	t = finish - start;
	cout << lanta << endl;
	cout << "迭代次数为：" << k << "   运行时间为：" << t << "ms" << endl << endl;
}
void exercise_6_1_3() {
	cout << "x^8+101x^7+208.01x^6+10891.01x^5+9802.08x^4+79108.9x^3-99902x^2+790x-1000=0的模最大根为：  " ;
	vector<vector<double>> A;
	vector<double> u(8, 0), A1 = { -101, -208.01 ,-10891.01 ,-9802.08 ,-79108.9 ,99902 ,-790 ,1000 };
	double lanta;
	int k;
	clock_t start, finish;
	double t;
	srand(time(0));
	u[0] = 1;
	A = make_zero_matrix(8, 8);
	A[0] = A1;
	for (int i = 1; i < 8; i++) A[i][i - 1] = 1;
	start = clock();
	lanta = Power_method(A, u, pow(10, -9),k);
	finish = clock();
	t = finish - start;
	cout << lanta << endl;
	cout << "迭代次数为：" << k << "   运行时间为：" << t << "ms" << endl << endl;
	//double x = -100;
	//cout << pow(x,8) + 101 * pow(x, 7) + 208.01 * pow(x, 6) + 10891.01 * pow(x, 5) + 9802.08 * pow(x, 4) + 79108.9 * pow(x, 3) - 99902 * pow(x, 2) + 790 * x - 1000;
}
void exercise_6_2_1() {
	vector<vector<double>> A;
	clock_t start, finish;
	double t;
	int k;
	srand(time(0));
	A = make_zero_matrix(41, 41);
	A[0][37] = A[0][40] = -1;
	for (int i = 1; i < 41; i++) A[i][i - 1] = 1;
	start = clock();
	A = Implicit_QR_algo(A, pow(10, -6), k);
	finish = clock();
	t = finish - start;
	cout_diag_characteristic_value(A);
	cout << endl;
	cout << "迭代次数为：" << k << "   运行时间为：" << t << "ms" << endl << endl;
}

void exercise_6_2_2() {
	vector<vector<double>> A={{9.1,3.0,2.6,4.0},{4.2,5.3,4.7,1.6},{3.2,1.7,9.4,0.8},{6.1,4.9,3.5,6.2}};
	clock_t start, finish;
	double t;
	int k;
	srand(time(0));
	//cout_Matrix(A);
	for (int i = 1; i < 4; i++) {
		A[2][3] += 0.1;
		start = clock();
		auto B = Implicit_QR_algo(A, pow(10, -6), k);
		finish = clock();
		t = finish - start;
		//cout_Matrix(A);
		cout_diag_characteristic_value(B);
		cout << endl;
		cout << "迭代次数为：" << k << "   运行时间为：" << t << "ms" << endl << endl;
	}
	
}


























