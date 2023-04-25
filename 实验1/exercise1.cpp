#include <iostream>
#include <vector>
#include"Function.h"
#include"exercise1.h"
#include<time.h>
using namespace std;
//ϰ��1.1
void exercise_1_1(void) {
	cout << "ϰ��1.1:"<<"\n";
	clock_t start, finish;
	double time;
	int m = 84;
	vector<vector<double>> A;
	vector<double> b(m, 7),x(m,1);
	
	A = make_zero_matrix(m, m);
	A[0][0] = 6;
	for (int i = 1; i < m; i++) {
		A[i][i] = 6;
		A[i - 1][i] = 1;
		A[i][i - 1] = 8;
		b[i] = 15;
	}
	b[m - 1] = 14;//��һ��ľ���A,����b����
	cout << "ȫ��Ԫ��ȥ���Ľ⣺" << "\n";
	start = clock();
	cout_Hvector(Gauss_all(A, m, b), m);
	finish = clock();
	time = finish - start;
	cout << "\n" << "��";
	cout<< Vector_norm((Vector_Subtracition(x, Gauss_all(A, m, b),m)),m,2) << "\n";
	cout << "the run time is :" << time << "ms" << "\n";

	cout << "����Ԫ��ȥ���Ľ⣺" << "\n";
	start = clock();
	cout_Hvector(Gauss_col(A, m, b), m);
	finish = clock();
	time = finish - start;
	cout << "\n" << "��";
	cout<<Vector_norm((Vector_Subtracition(x, Gauss_col(A, m, b),m)), m, 2) << "\n";
	cout << "the run time is :" << time << "ms" << "\n";

	cout << "��ѡ��Ԫ��ȥ���Ľ⣺" << "\n";
	start = clock();
	cout_Hvector(Gauss_no(A, m, b), m);
	finish = clock();
	time = finish - start;
	cout << "\n" << "��";
	cout << Vector_norm((Vector_Subtracition(x, Gauss_no(A, m, b),m)), m, 2) << "\n";
	cout << "the run time is :" << time << "ms" << "\n";
}

//ϰ��1.2.��1��+ϰ��1.3
void exercise_1_2_1(void) {
	cout << "ϰ��1.2.(1)+ϰ��1.3:" << "\n";
	clock_t start, finish;
	double t;
	int m = 100;
	vector<vector<double>> A,x,B,y;
	vector<double> b(m,0),q(m,0);
	srand(time(0));
	A = make_zero_matrix(m, m);
	y = make_zero_matrix(m, 1);
	x = make_zero_matrix(1, m);
	B = make_zero_matrix(m, 1);
	A[0][0] = 10;
	for (int i = 1; i < m; i++) {
		A[i][i] = 10;
		A[i - 1][i] = 1;
		A[i][i - 1] = 1;
	}
	for (int i = 0; i < m; i++) 
		y[i][0] = rand();
	B = Matrix_Mul(A, y, m, m, 1);
	for (int i = 0; i < m; i++) {
		b[i] = B[i][0];
		x[0][i] = y[i][0];
		q[i] = x[0][i];
	}                                //�ڶ���(1)�ľ���A,����b����,���������x��Ȼ��õ�b�������뾫ȷ��Ƚ�
	//cout_Matrix(A, m, m);
	//t_Lvector(b, m);
	cout << "��ʼ��b��" << "\n";
	cout_Hvector(b, m);
	cout << "\n";

	//cout << "��ȷ�⣺" << "\n";
	//cout_Hvector(q, m);
	//cout << "\n";

	cout << "ȫ��Ԫ��ȥ���Ľ⣺" << "\n";
	start = clock();
	cout_Hvector(Gauss_all(A, m, b), m);
	finish = clock();
	t = finish - start;
	cout << "��";
	cout << Vector_norm((Vector_Subtracition(q, Gauss_all(A, m, b),m)), m, 2)<< "\n";
	cout << "the run time is :" << t << "ms" << "\n";

	cout << "����Ԫ��ȥ���Ľ�" << "\n";
	start = clock();
	cout_Hvector(Gauss_col(A, m, b), m);
	finish = clock();
	t = finish - start;
	cout << "��";
	cout << Vector_norm((Vector_Subtracition(q, Gauss_col(A, m, b),m)), m, 2) << "\n";
	cout << "the run time is :" << t << "ms" << "\n";

	cout << "��ѡ��Ԫ��ȥ���Ľ�" << "\n";
	start = clock();
	cout_Hvector(Gauss_no(A, m, b), m);
	finish = clock();
	t = finish - start;
	cout << "��";
	cout << Vector_norm((Vector_Subtracition(q, Gauss_no(A, m, b),m)), m, 2) << "\n";
	cout << "the run time is :" << t << "ms" << "\n";

	cout << "ƽ�������Ľ�" << "\n";
	start = clock();
	cout_Hvector(Square_LLT(A, m, b), m);
	finish = clock();
	t = finish - start;
	cout << "��";
	cout << Vector_norm((Vector_Subtracition(q, Square_LLT(A, m, b),m)), m, 2) << "\n";
	cout << "the run time is :" << t << "ms" << "\n";

	cout << "�Ľ���ƽ�������Ľ�" << "\n";
	start = clock();
	cout_Hvector(Square_LDLT(A, m, b), m);
	finish = clock();
	t = finish - start;
	cout << "��";
	cout << Vector_norm((Vector_Subtracition(q, Square_LDLT(A, m, b),m)), m, 2) << "\n";
	cout << "the run time is :" << t << "ms" << "\n";
}

//ϰ��1.2.��2��+ϰ��1.3
void exercise_1_2_2(void){
	cout << "ϰ��1.2.(2)+ϰ��1.3:" << "\n";
	clock_t start, finish;
	double time;
	int m = 13;
	vector<vector<double>> A;
	vector<double> b(m),x(m,1);
	A = make_zero_matrix(m, m);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			A[i][j] =  1.0 / (i + j + 1);
			//A[i][j] =  1;
		}
	}
	for (int i = 0; i < m; i++) {
		b[i] = 0;
		for (int j = 0; j < m; j++) {
			b[i] +=  1.0 / (i + j + 1);
			//b[i] +=  1;
		}	
	}     //�ڶ���(2)�ľ���A,����b����
	//cout << "A" << "\n";
	//cout_Matrix(A, m, m);
	//cout_Hvector(b, m);
	cout << "ȫ��Ԫ��ȥ���Ľ⣺" << "\n";
	start = clock();
	cout_Hvector(Gauss_all(A, m, b), m);
	finish = clock();
	time = finish - start;
	cout << "\n";
	cout << "��";
	cout << Vector_norm((Vector_Subtracition(x, Gauss_all(A, m, b),m)), m, 2) << "\n";
	cout << "the run time is :" << time << "ms" << "\n";

	cout << "����Ԫ��ȥ���Ľ�" << "\n";
	start = clock();
	cout_Hvector(Gauss_col(A, m, b), m);
	finish = clock();
	time = finish - start;
	cout << "\n";
	cout << "��";
	cout << Vector_norm((Vector_Subtracition(x, Gauss_col(A, m, b),m)), m, 2) << "\n";
	cout << "the run time is :" << time << "ms" << "\n";

	cout << "��ѡ��Ԫ��ȥ���Ľ�" << "\n";
	start = clock();
	cout_Hvector(Gauss_no(A, m, b), m);
	finish = clock();
	time = finish - start;
	cout << "\n";
	cout << "��";
	cout << Vector_norm((Vector_Subtracition(x, Gauss_no(A, m, b),m)), m, 2) << "\n";
	cout << "the run time is :" << time << "ms" << "\n";

	cout << "ƽ�������Ľ�" << "\n";
	start = clock();
	cout_Hvector(Square_LLT(A, m, b), m);
	finish = clock();
	time = finish - start;
	cout << "\n";
	cout << "��";
	cout << Vector_norm((Vector_Subtracition(x, Square_LLT(A, m, b),m)), m, 2) << "\n";
	cout << "the run time is :" << time << "ms" << "\n";

	cout << "�Ľ���ƽ�������Ľ�" << "\n";
	start = clock();
	cout_Hvector(Square_LDLT(A, m, b), m);
	finish = clock();
	time = finish - start;
	cout << "\n";
	cout << "��";
	cout << Vector_norm((Vector_Subtracition(x, Square_LDLT(A, m, b),m)), m, 2) << "\n";
	cout << "the run time is :" << time << "ms" << "\n";
}
 








































