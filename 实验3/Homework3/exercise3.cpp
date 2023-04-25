#include <iostream>
#include <vector>
#include"Function.h"
#include"exercise3.h"
#include<time.h>
//#include<random>
#include<cstdlib>

using namespace std;
//习题3第一问
void exercise_3_1(void) {
	cout << "习题3.1:" << "\n";
	clock_t start, finish;
	double t;
	srand(time(0));
	cout << "QR分解求解第一章第一问方程组：" << endl;
	int m1 = 54;
	vector<vector<double>> A1;
	vector<double> b1(m1,7), d(m1,0),x1(m1,0),X1(m1,1);
	A1 = make_zero_matrix(m1, m1);
	A1[0][0] = 6;
	for (int i = 1; i < m1; i++) {
		A1[i][i] = 6;
		A1[i - 1][i] = 1;
		A1[i][i - 1] = 8;
		b1[i] = 15;
	}
	b1[m1 - 1] = 14;//第一章第一题的矩阵A,向量b输入
	start = clock();
	x1 = Q_R_solve(A1, b1, m1, m1);
	finish = clock();
	t = finish - start;
	cout <<"阶为" << m1 << "的解为：";
	cout_Hvector(x1, m1);
	cout <<"\n"<< "误差：";
	cout << Vector_norm((Vector_Subtracition(X1, x1, m1)), m1, 2) << "\n";
	cout <<  "the run time is : " << t << "ms" << "\n"<<endl;

	cout << "QR分解求解第一章第二问第一个方程组：" << endl;
	int m2 = 100;
	vector<vector<double>> A2;
	vector<double> b2(m2, 0),y(m2,0),x2(m2,0);
	A2 = make_zero_matrix(m2, m2);
	A2[0][0] = 10;
	for (int i = 1; i < m2; i++) {
		A2[i][i] = 10;
		A2[i - 1][i] = 1;
		A2[i][i - 1] = 1;
	}
	for (int i = 0; i < m2; i++)
		y[i] = rand();   
	b2 = Matrix_Mul_Vector(A2, y, m2, m2);//第一章第二题（1）的矩阵A,向量b输入
	/*cout_Matrix(A2, m2, m2);
	cout_Hvector(b2, m2);*/
	start = clock();
	x2 = Q_R_solve(A2, b2, m2, m2);
	finish = clock();
	t = finish - start;
	cout << "阶为" << m2 << "的解为：";
	cout_Hvector(x2, m2);
	cout <<"\n"<< "精确解为：";
	cout_Hvector(y, m2);
	cout << "\n" << "误差：";
	cout << Vector_norm((Vector_Subtracition(x2, y, m2)), m2, 2) ;
	cout << "\n" << "the run time is : " << t << "ms" << "\n"<<endl;

	cout << "QR分解求解第一章第二问第二个方程组：" << endl;
	int m3 = 40;
	vector<vector<double>> A3;
	vector<double> b3(m3), x3(m3, 1),X3(m3,1);
	A3 = make_zero_matrix(m3, m3);
	for (int i = 0; i < m3; i++) {
		for (int j = 0; j < m3; j++) {
			A3[i][j] = 1.0 / (i + j + 1);
			//A3[i][j] =  1;
		}
	}
	for (int i = 0; i < m3; i++) {
		b3[i] = 0;
		for (int j = 0; j < m3; j++) {
			b3[i] += 1.0 / (i + j + 1);
			//b3[i] +=  1;
		}
	}
	start = clock();
	x3 = Q_R_solve(A3, b3, m3, m3);
	finish = clock();
	t = finish - start;
	cout << "阶为" << m3 << "的Hilbert方程组的解为：";
	cout_Hvector(x3, m3);
	cout << "\n" << "误差：";
	cout << Vector_norm((Vector_Subtracition(X3, x3, m3)), m3, 2) ;
	cout << "\n" << "the run time is : " << t << "ms" << "\n" << endl;
}

//习题3第二问
void exercise_3_2(void) {
	cout << "习题三第二问："<< endl;
	vector<vector<double>> A;
	vector<double> y = {1.00,0.8125,0.75,1.00,1.3125,1.75,2.3125},x(3,0),X(3,1),t = {-1.00,-0.75,-0.5,0,0.25,0.5,0.75};
	double s;
	clock_t start, finish;
	srand(time(0));
	A = make_zero_matrix(7, 3);
	for (int i = 0; i < 7; i++) {
		A[i][0] = t[i] * t[i];
		A[i][1] = t[i];
		A[i][2] = 1;
	}
	//cout_Hvector(y, 7);
	//cout << endl;
	//cout_Matrix(A, 7, 3);
	start = clock();
	x = Q_R_solve(A, y, 7, 3);
	finish = clock();
	s = finish - start;
	cout  <<"解为：";
	cout_Hvector(x, 3);
	cout << "\n" << "误差：";
	cout << Vector_norm((Vector_Subtracition(x, X, 3)), 3, 2) << "\n";
	cout << "the run time is : " << s << "ms" << "\n" << endl;
}


//习题3第三问
void exercise_3_3(void) {
	cout << "习题三第三问：" << endl;
	vector<vector<double>> A =
	{ {1,4.9176, 1, 3.472, 0.998, 1, 7, 4, 42, 3, 1, 0},
	{1,5.0208, 1, 3.531, 1.5, 2, 7, 4, 62, 1, 1, 0},
	{1,4.5429, 1, 2.275, 1.175, 1, 6, 3, 40,  2, 1, 0},
	{1,4.5573, 1, 4.05, 1.232, 1, 6, 3, 54, 4, 1, 0},
	{1,5.0597, 1, 4.455, 1.121, 1, 6, 3, 42, 3, 1, 0},
	{1,3.891, 1, 4.455, 0.988, 1, 6, 3, 56, 2, 1, 0},
	{1,5.898, 1, 5.85, 1.24, 1, 7, 3, 51, 2, 1,  1},
	{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
	{1,15.4202, 2.5,  9.8, 3.42, 2, 10, 5, 42, 2, 1, 1},
	{1,14.4598, 2.5, 12.8, 3, 2, 9, 5, 14, 4, 1, 1},
	{1,5.8282, 1, 6.435, 1.225, 2, 6, 3, 32, 1, 1, 0},
	{1,5.3003, 1, 4.9883, 1.552, 1, 6, 3, 30, 1, 2, 0},
	{1,6.2712, 1, 5.52, 0.975, 1, 5, 2, 30, 1, 2, 0},
	{1,5.9592, 1, 6.666, 1.121, 2, 6, 3, 32, 2, 1, 0},
	{1,5.05, 1, 5, 1.02, 0, 5, 2, 46, 4, 1, 1},
	{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
	{1,8.2464, 1.5, 5.15, 1.664, 2, 8, 4, 50, 4, 1, 0},
	{1,6.6969, 1.5, 6.092, 1.488, 1.5, 7, 3, 22, 1, 1, 1},
	{1,7.7841, 1.5, 7.102, 1.376, 1, 6, 3, 17, 2, 1, 0},
	{1,9.0384, 1, 7.8, 1.5, 1.5, 7, 3, 23, 3, 3, 0},
	{1,5.9894, 1, 5.52, 1.256, 2, 6, 3, 40, 4, 1, 1},
	{1,7.5422, 1.5, 4, 1.69, 1, 6, 3, 22, 1, 1, 0},
	{1,8.7951, 1.5, 9.89, 1.82, 2, 8, 4, 50, 1, 1, 1},
	{1,6.0931, 1.5, 6.7265, 1.652, 1, 6, 3, 44, 4, 1, 0},
	{1,8.3607, 1.5, 9.15, 1.777, 2., 8, 4, 48, 1, 1, 1},
	{1,8.14, 1, 8, 1.504, 2, 7, 3, 3, 1, 3, 0},
	{1,9.1416, 1.5, 7.3262, 1.831, 1.5, 8, 4, 31, 4, 1, 0},
	{1,12, 1.5, 5, 1.2, 2, 6, 3, 30, 3, 1, 1} };
	vector<double> b =
	{ 25.9, 29.5, 27.9, 25.9, 29.9, 29.9, 30.9,
	28.9, 84.9, 82.9, 35.9, 31.5, 31.0, 30.9,
	30.0, 28.9, 36.9, 41.9, 40.5, 43.9, 37.5,
	37.9, 44.5, 37.9, 38.9, 36.9, 45.8, 41.0 },x(12,0);
	double s;
	clock_t start, finish;
	srand(time(0));
	start = clock();
	x = Q_R_solve(A, b, 28, 12);
	finish = clock();
	s = finish - start;
	cout << "解为：";
	cout_Hvector(x, 12);
	//cout_Hvector(Matrix_Mul_Vector(A, x, 28, 12), 28);
	cout << "\n" << "误差：";
	cout << Vector_norm((Vector_Subtracition(Matrix_Mul_Vector(A,x,28,12), b, 28)), 28, 2) << "\n";
	cout << "the run time is : " << s << "ms" << "\n" << endl;
}









































