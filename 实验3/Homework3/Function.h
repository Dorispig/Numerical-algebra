#pragma once
#include<iostream>
#include <vector>
#include<math.h>
#define double long double
using namespace std;

vector<vector<double>> Matrix_add(vector<vector<double>> a, vector<vector<double>> b);                 //矩阵相加
vector<vector<double>> Matrix_Mul(vector<vector<double>> a, vector<vector<double>> b,const int& m, const int& n, const int& p); //矩阵相乘(默认可以相乘，主函数里判断）
vector<double> Matrix_Mul_Vector(vector<vector<double>> A, vector<double> x, int m, int n);            //矩阵乘向量
vector<double> Vector_Mul_Matrix(vector<double> x ,vector<vector<double>> A, int m, int n);            //向量乘矩阵
vector<vector<double>> LVector_Mul_HVector(vector<double> x, vector<double> y, int m, int n);    //列向量乘行向量

vector<vector<double>> Matrix_T(vector<vector<double>> a,const int& m, const int& n);                  //求矩阵转置
vector<vector<double>> Matrix_Mul_num(vector<vector<double>> a, double b);                             //矩阵数乘
vector<vector<double>> Matrix_Divid_num(vector<vector<double>> a, double b);                           //矩阵除以一个数
vector<vector<double>> Matrix_Subtracition(vector<vector<double>> a, vector<vector<double>> b,int m,int n);        //矩阵相减,a-b
vector<double> Vector_Subtracition(vector<double> a, vector<double> b,int m);                          //向量相减
vector<vector<double>> make_zero_matrix(int a, int b);                                                 //创建0矩阵
vector<vector<double>> I(int b);                                                                         //创建b阶单位矩阵
double det(vector<vector<double>> A,int n);                                                            //化为上三角矩阵求det（A）
vector<vector<double>> read_Matrix(int m, int n);                                                      // 读入m * n矩阵（）
vector<vector<double>> Matrix_block(vector<vector<double>> A, int start_x, int start_y, int block_x, int block_y);//分块矩阵(从start_x,start_y开始，大小为block_x * block_y的分块矩阵)
vector<double> Upper_triangle(int n,vector<vector<double>> U, vector<double> b);                       //求解上三角阵方程组Ux = b
vector<double> Lower_triangle(int n,vector<vector<double>> L, vector<double> b);                       //求解下三角阵方程组Lx = b
void cout_Matrix(vector<vector<double>> A,int m,int n);                                                //输出矩阵
void cout_Lvector(vector<double> b, int n);                                                            //输出列向量
void cout_Hvector(vector<double> b, int n);                                                            //输出行向量
void swap_row(vector<vector<double>>& A, int m, int n, int k, int p);                                  //交换m * n阶方阵第k行和第p行
void swap_col(vector<vector<double>>& A, int m, int n, int k, int p);                                  //交换m * n阶方阵第k列和第p列
void swap_vector(vector<double>& b, int k, int p);                                                     //交换向量b的第k,p个元素
void max_subMatrix(vector<vector<double>> A, int j, int t, int m, int n, int& p, int& q);              //找某个子矩阵（ij--mn）里绝对值的最大元素位置
void getLU(vector<vector<double>> A, const int& n, vector<vector<double>>& L, vector<vector<double>>& U);//取L,U矩阵
vector<double> Gauss_no(vector<vector<double>> A,const int &n,vector<double> b);                       //不选主元的Gauss消去法
vector<double> Gauss_all(vector<vector<double>> A, const int& n, vector<double> b);                    //全主元Gauss消去法
vector<double> Gauss_col(vector<vector<double>> A, const int& n,vector<double> b);                     //列主元Gauss消去法
//void Gauss_col_LU(vector<vector<double>> A, const int& n, vector<double>& L, vector<double>& U);       //列主元三角分解求L，U
vector<double> Square_LLT(vector<vector<double>> A,const int& n, vector<double> b);                    //平方根消去法
vector<double> Square_LDLT(vector<vector<double>> A, const int& n, vector<double> b);                  //改进的平方根消去法
double Vector_norm( vector<double> b, int n,double k);                                                 //向量范数
double Matrix_norm(vector<vector<double>> A, int m, int n,int k);                                      //矩阵范数
double Matrix_norm_1_pro(vector<vector<double>> A,int n);                                              //优化法估计矩阵A的逆的无穷范数
vector<vector<double>> sign_Matrix(vector<vector<double>> w, int m, int n);                            //判断矩阵元素符号
double K(vector<vector<double>> A, int n,int k);                                                       //求A的条件数  
void House(vector<double> x, vector<double> &v, double &beta, int n);                                  //计算Householder变换
void givens(double c, double s, double a, double b);                                                   //计算givens变换
vector<vector<double>> Q_R(vector<vector<double>> A, vector<double>& d, int m, int n);                         //计算A的QR分解
vector<double> Q_R_solve(vector<vector<double>> A, vector<double> b, int m, int n);                      //用QR分解解方程




