#include <iostream>
#include <vector>
#include"Function.h"
#include"exercise4.h"
#include<time.h>
//#include<random>
#include<cstdlib>

using namespace std;


double f(double x,double e,double a) {
	return (1.0 - a) / (1.0 - exp( -1 / e)) * (1 - exp( -x / e)) + a * x;
}


void exercise_4_1(double e, double a, int n,double w) {
	n = n - 1;
	vector<vector<double>> A;
	vector<double> b(n), x(n),x1(n,1),x2(n,0);
	double  h = 1.0 / (n+1);
	clock_t start, finish;
	double t;
	srand(time(0));
	A = make_zero_matrix(n, n);
	A[0][0] = -2 * e - h;
	b[0] = a * h * h;
	x[0] = f(h, e, a); 
	for (int i = 1; i < n; i++) {
		A[i][i] = A[0][0];
		A[i][i - 1] = e;
		A[i - 1][i] = e + h;
		b[i] = b[0];
		x[i] = f((i+1) * h, e, a);
	}
	b[n - 1] -= (e + h);
	cout << "e = " << e << " " << "a = " << a << " " << "n = " << n+1 << "时三种方法的结果如下：" << endl;
	cout << "\n";
	//cout_Hvector(x); cout << endl;

	                      //Jacobi迭代法
	int Sum_Jacobi = 0;
	cout << "Jacobi法结果：" << endl;
	//cout << Vector_norm(Vector_Subtracition(x1, x2), -1) << endl;
	start = clock();
	while (Vector_norm(Vector_Subtracition(x1, x2), -1) > pow(10,-6)){
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
	cout << "与精确解的误差为：" << Vector_norm(Vector_Subtracition(x1, x), -1) << endl;

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
	cout << "与精确解的误差为：" << Vector_norm(Vector_Subtracition(x1, x), -1) << endl;

	                     //SOR迭代法
	for (int i = 0; i < n; i++) {
		x1[i] = 1;
		x2[i] = 0;
	}
	cout << "\nSOR迭代法结果：" << endl;
	int Sum_SOR = 0;
	start = clock();
	while (Vector_norm(Vector_Subtracition(x1, x2), -1) > pow(10, -6)) {
		x1 = x2;
		x2 = SOR_Iterative_method(A, b, x2,w);
		//cout_Hvector(x2); cout << endl;
		Sum_SOR++;
	} 
	finish = clock();
	t = finish - start;
	cout << "解为：";
	cout_Hvector(x2); cout << endl;
	cout << "松弛因子为：" << w << endl;
	cout << "迭代次数为：" << Sum_SOR << endl;
	cout << "运行时间为：" << t << "ms" << endl;
	cout << "与精确解的误差为：" << Vector_norm(Vector_Subtracition(x1, x), -1) << endl;






}


void exercise_4_2(int N,double w) {
	vector<vector<double>> F, G, U1, U2, U3,U;
	double h = 1.0 / N;
	clock_t start, finish;
	double t;
	srand(time(0));
	F = make_zero_matrix(N + 1, N + 1);
	G = make_zero_matrix(N + 1, N + 1);
	U1 = make_zero_matrix(N + 1, N + 1);
	U2 = make_zero_matrix(N + 1, N + 1);
	U3 = make_zero_matrix(N + 1, N + 1);
	U = make_zero_matrix(N + 1, N + 1);
	for (int i = 0; i <= N; i++) {
		for (int j = 0; j <= N; j++) {
			F[i][j] = h * h * (i*h+j*h);
			G[i][j] = h * h * exp(i * h * j * h);
			if ((i == 0) || (i == N) || (j == 0) || (j == N)) {
				U1[i][j] = 1;
				U2[i][j] = 1;
				U3[i][j] = 1;
			}
			//else {
			//	U1[i][j] = 0.5;
			//	U2[i][j] = 0.5;
			//	U3[i][j] = 0.5;
			//}
		}
	}
	//Jacobi迭代法
	int Sum_Jacobi = 1;
	cout << "Jacobi法结果：" << endl;
	//cout << Vector_norm(Vector_Subtracition(x1, x2), -1) << endl;
	start = clock();
	U = Jacobi_Iterative_method2(U1, F, G);
	while ((Matrix_norm(Matrix_Subtracition(U,U1),-1)) > pow(10, -7)) {
		U = U1;
		U1 = Jacobi_Iterative_method2(U1, F, G);
		//cout_Hvector(x2); cout << endl;
		Sum_Jacobi++;
		//cout << Vector_norm(Vector_Subtracition(x1, x2), -1) << endl;
	}
	finish = clock();
	t = finish - start;
	cout << "解为：" << endl;
	cout_Matrix(Matrix_block(U, 1, 1, N - 1, N - 1));
	cout << "迭代次数为：" << Sum_Jacobi << endl;
	cout << "运行时间为：" << t << "ms" << endl;

	//G_S迭代法
	int G_S_Jacobi = 1;
	cout << "G-S迭代法结果：" << endl;
	//cout << Vector_norm(Vector_Subtracition(x1, x2), -1) << endl;
	start = clock();
	U = G_S_Iterative_method2(U2, F, G);
	while ((Matrix_norm(Matrix_Subtracition(U, U2), -1)) > pow(10, -7)) {
		U = U2;
		U2 = G_S_Iterative_method2(U2,F,G);
		//cout_Hvector(x2); cout << endl;
		G_S_Jacobi++;
		//cout << Vector_norm(Vector_Subtracition(x1, x2), -1) << endl;
	}
	finish = clock();
	t = finish - start;
	cout << "解为：" << endl;
	cout_Matrix(Matrix_block(U, 1, 1, N - 1, N - 1));
	cout << "迭代次数为：" << G_S_Jacobi << endl;
	cout << "运行时间为：" << t << "ms" << endl;

	//SOR迭代法
	int SOR_Jacobi = 1;
	cout << "SOR迭代法结果：" << endl;
	//cout << Vector_norm(Vector_Subtracition(x1, x2), -1) << endl;
	start = clock();
	U = SOR_Iterative_method2(U3, F, G,w);
	while ((Matrix_norm(Matrix_Subtracition(U, U3), -1)) > pow(10, -7)) {
		U = U3;
		U3 = SOR_Iterative_method2(U3, F, G,w);
		//cout_Hvector(x2); cout << endl;
		SOR_Jacobi++;
		//cout << Vector_norm(Vector_Subtracition(x1, x2), -1) << endl;
	}
	finish = clock();
	t = finish - start;
	cout << "解为：" << endl;
	cout_Matrix(Matrix_block(U,1,1,N-1,N-1));
	cout << "松弛因子为：" << w << endl;
	cout << "迭代次数为：" << SOR_Jacobi << endl;
	cout << "运行时间为：" << t << "ms" << endl;
}







































