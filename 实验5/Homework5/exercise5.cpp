#include <iostream>
#include <vector>
#include"Function.h"
#include"exercise5.h"
#include<time.h>
//#include<random>
#include<cstdlib>

using namespace std;


void exercise_5_1(double W) {
	int n = 20, N = (n - 1) * (n - 1), k = 1, n0 = n - 1,P,Q;
	double a,beta,h =1.0/n,w;
	vector<vector<double>> A;
	vector<double>x(N,0.5), x1(N), x2(N), r(N), p(N),b(N);
	A = make_zero_matrix(N, N);
	clock_t start, finish;
	double t;
	srand(time(0));
	//初始化A,b
	w = h * h / 4;
	A[0][0] = 1 + w;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A[i][i] = A[0][0];
			if (abs(i - j) == n0) A[i][j] = -1 / 4.0;
			if (i < N - 1) {
				if ((i + 1) % n0 != 0) {
					A[i][i + 1] = -0.25;
					A[i + 1][i] = -0.25;
				}
			}
		}
	}//初始化A
	for (int i = 0; i < N; i++){
		if (i == 0) {
			b[0] = w * sin(h * h) + 0.5 * 4 * w;
			b[n0 - 1] = w * sin(n0 * h * h) + 0.25 * (1 + h * h + N * h * h);
			for (int j = 1; j < n0 - 1; j++) b[j] = w * sin((j+1) * h * h) + 0.25 * (j+1) * h * (j+1) * h;
			i = n0-1;
		}
		else if (i == n0*(n0 - 1)) {
			b[n0 * (n0 - 1)] = w * sin(h * n0 * h) + 0.25 * (1 + h * h + N * h * h);
			b[n0*n0-1] = w * sin(N * h * h) + 0.25 * (2 + 2 * N * h * h);
			for (int j = 1; j < n0 - 1; j++) b[n0 * (n0 - 1) + j] = w * sin((j + 1) * n0 * h * h) + 0.25 * (1 + (j + 1) * (j + 1) * h * h);
			i = N;
		}
		else {
			P = i % n0 + 1; Q = i / n0 + 1;
			if (P == 1) b[i] = w * sin(P * Q * h * h) + 0.25 * Q * Q * h * h;
			else if (P == n0) b[i] = w * sin(P * Q * h * h) + 0.25 * (1 + Q * Q * h * h);
			else b[i] = w * sin(P * Q * h * h);
		}
	}//初始化b;
	//cout << w * sin(n0 * h * h) + 0.25 * (1 + h * h + N * h * h) << endl;
	//cout_Matrix(A);
	//cout_Hvector(b); cout << endl;
	//r = Vector_Subtracition(b, Matrix_Mul_Vector(A, x));                            //r0 = b-Ax0
	//p = r;                                                                          //p0 = r0
	//a = HVector_Mul_LVector(r, r) / HVector_Mul_LVector(p, Matrix_Mul_Vector(A, p));//a0 = r0Tr0/p0TAp0
	//x1 = Vector_Add(x, Vector_Mul_num(p,a));                                        //x1 = x0 + a0p0
	//r = Vector_Subtracition(b, Matrix_Mul_Vector(A, x1));                           //r1 = b-Ax1
	//beta =  - HVector_Mul_LVector(r, Matrix_Mul_Vector(A, p)) / HVector_Mul_LVector(p, Matrix_Mul_Vector(A, p));//beta0=r1TAp0/p0TAp0
	//p = Vector_Add(r, Vector_Mul_num(p, beta));                                        //p1 = r1 +beta0*p0
	//x2 = x;
	start = clock();
	//while ((Vector_norm(Vector_Subtracition(x1, x2), -1)) > pow(10,-7)) {
	//	x2 = x1;
	//	conjugate_grad_method_1(A, b, x1, p, r);
	x = conjugate_grad_method(A, b, pow(10, -7), k);
	//	k++;
	//}
	finish = clock();
	t = finish - start;
	cout << "共轭梯度法的解为：";
	cout_Hvector(x); cout << endl;
	cout << "迭代次数为：" << k << endl;
	cout << "运行时间为：" << t << "ms" << endl;

	//SOR迭代法
	for (int i = 0; i < n; i++) {
		x1[i] = 1;
		x2[i] = 0;
	}
	cout << "\nSOR迭代法结果：" << endl;
	int Sum_SOR = 0;
	start = clock();
	while (Vector_norm(Vector_Subtracition(x1, x2), -1) > pow(10, -7)) {
		x1 = x2;
		x2 = SOR_Iterative_method(A, b, x2, W);
		//cout_Hvector(x2); cout << endl;
		Sum_SOR++;
	}
	finish = clock();
	t = finish - start;
	cout << "解为：";
	cout_Hvector(x2); cout << endl;
	cout << "松弛因子为：" << W << endl;
	cout << "迭代次数为：" << Sum_SOR << endl;
	cout << "运行时间为：" << t << "ms" << endl;

	
}

void exercise_5_2(int n) {
	vector<vector<double>> A;
	vector<double>x(n, 0.5), x1(n), x2(n), r(n), p(n), b(n),X(n,1.0/3);
	clock_t start, finish;
	double t, a, beta;
	int k = 0;
	srand(time(0));
	A = make_zero_matrix(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) A[i][j] = 1.0 / (i + j + 1);
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			b[i] += A[i][j]/3.0;
	}
	//cout_Matrix(A);
	//cout_Hvector(b);
	//r = Vector_Subtracition(b, Matrix_Mul_Vector(A, x));                            //r0 = b-Ax0
	//p = r;                                                                          //p0 = r0
	//a = HVector_Mul_LVector(r, r) / HVector_Mul_LVector(p, Matrix_Mul_Vector(A, p));//a0 = r0Tr0/p0TAp0
	//x1 = Vector_Add(x, Vector_Mul_num(p, a));                                        //x1 = x0 + a0p0
	//r = Vector_Subtracition(b, Matrix_Mul_Vector(A, x1));                           //r1 = b-Ax1
	//beta = -HVector_Mul_LVector(r, Matrix_Mul_Vector(A, p)) / HVector_Mul_LVector(p, Matrix_Mul_Vector(A, p));//beta0=r1TAp0/p0TAp0
	//p = Vector_Add(r, Vector_Mul_num(p, beta));                                        //p1 = r1 +beta0*p0
	//x2 = x;
	start = clock();
	//while ((Vector_norm(Vector_Subtracition(x1, x2), -1)) > pow(10, -7)) {
		//x2 = x1;
		//conjugate_grad_method_1(A, b, x1, p, r);
		//k++;
		x = conjugate_grad_method(A, b, pow(10, -7), k);
	//}
	finish = clock();
	t = finish - start;
	cout << "阶为" << n << "的Hilibert矩阵的共轭梯度法的解为：";
	cout_Hvector(x); cout << endl;
	cout << "迭代次数为：" << k << endl;
	cout << "运行时间为：" << t << "ms" << endl;
	cout << "与精确解的误差为：" << Vector_norm(Vector_Subtracition(X, x), -1) << endl << endl;
}

void exercise_5_3(void) {
	//int n = 5;
	//vector<vector<double>> A = { {10,1,2,3,4},{1,9,-1,2,-3},{2,-1,7,3,-5},{3,2,3,12,-1},{4,-3,-5,-1,15}};
	//vector<double>x(n, 0.5), x1(n), x2(n), r(n), p(n), b={12,-27,14,-17,12}, X(n, 1);
	//clock_t start, finish;
	//double t, a, beta;
	//int k = 0;
	//srand(time(0));
	//r = Vector_Subtracition(b, Matrix_Mul_Vector(A, x));                            //r0 = b-Ax0
	//p = r;                                                                          //p0 = r0
	//a = HVector_Mul_LVector(r, r) / HVector_Mul_LVector(p, Matrix_Mul_Vector(A, p));//a0 = r0Tr0/p0TAp0
	//x1 = Vector_Add(x, Vector_Mul_num(p, a));                                        //x1 = x0 + a0p0
	//r = Vector_Subtracition(b, Matrix_Mul_Vector(A, x1));                           //r1 = b-Ax1
	//beta = -HVector_Mul_LVector(r, Matrix_Mul_Vector(A, p)) / HVector_Mul_LVector(p, Matrix_Mul_Vector(A, p));//beta0=r1TAp0/p0TAp0
	//p = Vector_Add(r, Vector_Mul_num(p, beta));                                        //p1 = r1 +beta0*p0
	//x2 = x;
	//start = clock();
	//while ((Vector_norm(Vector_Subtracition(x1, x2), -1)) > pow(10, -7)) {
	//	x2 = x1;
	//	conjugate_grad_method_1(A, b, x1, p, r);
	//	k++;
	//}
	//finish = clock();
	//t = finish - start;
	//cout << "5.3的解为：" << endl;
	//cout_Hvector(x1); cout << endl;
	//cout << "迭代次数为：" << k << endl;
	//cout << "运行时间为：" << t << "ms" << endl;
	//cout << "与精确解的误差为：" << Vector_norm(Vector_Subtracition(X, x1), -1) << endl;

	int n = 5;
	vector<vector<double>> A = { {10,1,2,3,4},{1,9,-1,2,-3},{2,-1,7,3,-5},{3,2,3,12,-1},{4,-3,-5,-1,15} };
	vector<double> b = { 12,-27,14,-17,12 },x(n),x1(n, 1), x2(n, 0);
	clock_t start, finish;
	double t, a, beta;
	int k = 0;
	srand(time(0));
	
	start = clock();
	x = conjugate_grad_method(A, b, pow(10, -7),k);
	finish = clock();
	t = finish - start;
	cout << "5.3的解为：" << endl;
	cout_Hvector(x); cout << endl;
	cout << "迭代次数为：" << k << endl;
	cout << "运行时间为：" << t << "ms" << endl << endl;

	//Jacobi迭代法
	int Sum_Jacobi = 0;
	cout << "Jacobi法结果：" << endl;
	//cout << Vector_norm(Vector_Subtracition(x1, x2), -1) << endl;
	start = clock();
	while (Vector_norm(Vector_Subtracition(x1, x2), -1) > pow(10, -6)) {
		x1 = x2;
		x2 = Jacobi_Iterative_method(A, b, x2);
		//cout_Hvector(x2); cout << endl;
		Sum_Jacobi++;
		//cout << Vector_norm(Vector_Subtracition(x1, x2), -1) << endl;
	}
	finish = clock();
	t = finish - start;
	cout << "解为：";
	cout_Hvector(x2); cout << endl;
	cout << "迭代次数为：" << Sum_Jacobi << endl;
	cout << "运行时间为：" << t << "ms" << endl;
	//cout << "与精确解的误差为：" << Vector_norm(Vector_Subtracition(x1, x), -1) << endl;

	//G-S迭代法
	for (int i = 0; i < n; i++) {
		x1[i] = 1;
		x2[i] = 0;
	}
	cout << "\nG-S迭代法结果：" << endl;
	int Sum_G_S = 0;
	start = clock();
	while (Vector_norm(Vector_Subtracition(x2, x1), -1) > pow(10, -6)) {
		x1 = x2;
		x2 = G_S_Iterative_method(A, b, x2);
		Sum_G_S++;
	}
	finish = clock();
	t = finish - start;
	cout << "解为：";
	cout_Hvector(x2); cout << endl;
	cout << "迭代次数为：" << Sum_G_S << endl;
	cout << "运行时间为：" << t << "ms" << endl;
	////cout << "与精确解的误差为：" << Vector_norm(Vector_Subtracition(x1, x), -1) << endl;
}































