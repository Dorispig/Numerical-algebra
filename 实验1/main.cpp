// 数值代数.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
#include <iostream>
#include <vector>
#include"Function.h"
#include"exercise1.h"
using namespace std;
//vector<vector<double>>
//上级习题1，求解84阶方程组

//int main()
//{
//    int m = 2;
//    vector<vector<double>> A;
//    vector<double> b(m);
////    vector<vector<double>>L;
////    vector<vector<double>>U;
//    A = make_zero_matrix(m, m);
////    L = make_zero_matrix(m, m);
////    U = make_zero_matrix(m, m);
////    A[0][0] = 6;
//    A[0][0] = 375; 
//    A[0][1] = 374; 
//    A[1][0] = 752 ;
//    A[1][1] = 750;
//    b[0] = 749.01;
//    b[1] = 1502.01;
////    /*for (int i = 1; i < m; i++) {
////        A[i][i] = 6;
////        A[i - 1][i] = 1;
////        A[i][i - 1] = 8;
////        b[i] = 15;
////    }
////    b[m - 1] = 14;*///第一题的矩阵A,向量b输入
////    //A[0][0] = 10;
////    //for (int i = 1; i < m; i++) {
////    //    A[i][i] = 10;
////    //    A[i - 1][i] = 1;
////    //    A[i][i - 1] = 1;
////    //    b[i] = 12;
////    //}
////    //b[m - 1] = 11;//第二题(1)的矩阵A,特殊的向量b输入
////    //
////    //cout_Lvector(b, m);
////    //Gauss_no(A, m, b);//不选主元的gauss消去法
////    //swap_row(A, m, m, 1, 2);
// //   b=Gauss_all(A, m, b);//全主元的gauss消去法
//    b = Gauss_col(A, m, b);//列主元的gauss消去法
////    Square_LLT(A, m, b);//平方根法
////    //Square_LDLT(A, m, b);//改进的平方根法
////    //getLU(A, m, L, U);
////    //Lower_triangle(m, L, b);
////    //Upper_triangle(m, U, b);
////    //cout_Matrix(L, m, m);
////    //cout_Matrix(U, m, m);
//    cout_Lvector(b, m);
////    //cout_Matrix(A,m,m);
////    //cout << p << q;
//    return 0;
//}
////    生成（a，b）内的随机实数：
////    #include<random>(此电脑好像不能用？)
////    #include<ctime>;
////    uniform_real_distribution<double> u(a,b);
////    default_random_engine e(time(NULL));
////    for(int i = 0;i<10;i++){
////        a[i] = u(e);
////    }
//
//
//
//
//
//
//
//
//
//
//
//// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
