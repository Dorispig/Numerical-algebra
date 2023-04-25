#pragma once
#include<iostream>
#include <vector>
#include<math.h>
#define double long double
using namespace std;

vector<vector<double>> Matrix_add(vector<vector<double>> a, vector<vector<double>> b);                 //矩阵相加
vector<vector<double>> Matrix_Mul(vector<vector<double>> a, vector<vector<double>> b);                 //矩阵相乘(默认可以相乘，主函数里判断）
vector<double> Matrix_Mul_Vector(vector<vector<double>> A, vector<double> x);                          //矩阵乘向量
vector<double> Vector_Mul_Matrix(vector<double> x ,vector<vector<double>> A);                          //向量乘矩阵
vector<vector<double>> LVector_Mul_HVector(vector<double> x, vector<double> y);                        //列向量乘行向量
double HVector_Mul_LVector(vector<double> x, vector<double> y);                                        //行向量乘列向量
vector<vector<double>> Matrix_T(vector<vector<double>> a);                                             //求矩阵转置
vector<vector<double>> Matrix_Mul_num(vector<vector<double>> a, double b);                             //矩阵数乘
vector<double> Vector_Mul_num(vector<double> a, double b);                                             //向量数乘
vector<vector<double>> Matrix_Divid_num(vector<vector<double>> a, double b);                           //矩阵除以一个数
vector<vector<double>> Matrix_Subtracition(vector<vector<double>> a, vector<vector<double>> b);        //矩阵相减,a-b
vector<double> Vector_Subtracition(vector<double> a, vector<double> b);                                //向量相减
vector<double> Vector_Add(vector<double> a, vector<double> b);                                         //向量相加
vector<vector<double>> make_zero_matrix(int a, int b);                                                 //创建0矩阵
vector<vector<double>> I_matrix(int b);                                                                //创建b阶单位矩阵
double det(vector<vector<double>> A);                                                                  //化为上三角矩阵求det（A）
vector<vector<double>> read_Matrix(int m, int n);                                                      // 读入m * n矩阵（）
vector<vector<double>> Matrix_block(vector<vector<double>> A, int start_x, int start_y, int block_x, int block_y);//分块矩阵(从start_x,start_y开始，大小为block_x * block_y的分块矩阵)
vector<vector<double>> Matrix_block_mul(vector<vector<double>> A, int start_x, int start_y, int block_x, int block_y,vector<vector<double>>B,int k);//将A中分块矩阵更新为乘B之后的值,k=-1左乘k=1右乘
vector<double> Upper_triangle(vector<vector<double>> U, vector<double> b);                             //求解上三角阵方程组Ux = b
vector<double> Lower_triangle(vector<vector<double>> L, vector<double> b);                             //求解下三角阵方程组Lx = b
void cout_Matrix(vector<vector<double>> A);                                                            //输出矩阵
void cout_Lvector(vector<double> b);                                                                   //输出列向量
void cout_Hvector(vector<double> b);                                                                   //输出行向量
void swap_row(vector<vector<double>>& A, int k, int p);                                                //交换m * n阶方阵第k行和第p行
void swap_col(vector<vector<double>>& A, int k, int p);                                                //交换m * n阶方阵第k列和第p列
void swap_vector(vector<double>& b, int k, int p);                                                     //交换向量b的第k,p个元素
void max_subMatrix(vector<vector<double>> A, int j, int t, int m, int n, int& p, int& q);              //找某个子矩阵（ij--mn）里绝对值的最大元素位置
void getLU(vector<vector<double>> A, vector<vector<double>>& L, vector<vector<double>>& U);            //取L,U矩阵
vector<double> Gauss_no(vector<vector<double>> A,vector<double> b);                                    //不选主元的Gauss消去法
vector<double> Gauss_all(vector<vector<double>> A,  vector<double> b);                                 //全主元Gauss消去法
vector<double> Gauss_col(vector<vector<double>> A, vector<double> b);                                  //列主元Gauss消去法
//void Gauss_col_LU(vector<vector<double>> A, const int& n, vector<double>& L, vector<double>& U);     //列主元三角分解求L，U
vector<double> Square_LLT(vector<vector<double>> A, vector<double> b);                                 //平方根消去法
vector<double> Square_LDLT(vector<vector<double>> A,  vector<double> b);                               //改进的平方根消去法
double Vector_norm( vector<double> b,double k);                                                        //向量范数
double Matrix_norm(vector<vector<double>> A,int k);                                                    //矩阵范数
double Matrix_norm_1_pro(vector<vector<double>> A);                                                    //优化法估计矩阵A的逆的无穷范数
vector<vector<double>> sign_Matrix(vector<vector<double>> w);                                          //判断矩阵元素符号
double K(vector<vector<double>> A,int k);                                                              //求A的条件数  
void House(vector<double> x, vector<double> &v, double &beta);                                         //计算Householder变换
void givens(double c, double s, double a, double b);                                                   //计算givens变换
vector<vector<double>> Q_R(vector<vector<double>> A, vector<double>& d);                               //计算A的QR分解
vector<double> Q_R_solve(vector<vector<double>> A, vector<double> b);                                  //用QR分解解方程
vector<double> Jacobi_Iterative_method(vector<vector<double>> A, vector<double> b, vector<double> x);  //Jacobi迭代法解方程
vector<double> G_S_Iterative_method(vector<vector<double>> A, vector<double> b, vector<double> x);     //G-S迭代法解方程
vector<double> SOR_Iterative_method(vector<vector<double>> A, vector<double> b, vector<double> x, double w);     //SOR迭代法解方程
vector<vector<double>> Jacobi_Iterative_method2(vector<vector<double>> U, vector<vector<double>> F, vector<vector<double>> G);          //Jacobi迭代法解方程，二维偏微分方程
vector<vector<double>> G_S_Iterative_method2(vector<vector<double>> A, vector<vector<double>> F, vector<vector<double>> G);               //G-S迭代法解方程,二维偏微分方程
vector<vector<double>> SOR_Iterative_method2(vector<vector<double>> A, vector<vector<double>> F, vector<vector<double>> G, double w);   //SOR迭代法解方程,二维偏微分方程
void conjugate_grad_method_1(vector<vector<double>> A, vector<double> b, vector<double> &x, vector<double> &p, vector<double> &r);       //共轭梯度法
vector<double> conjugate_grad_method(vector<vector<double>> A, vector<double> b,double w, int &k);//共轭梯度法求解,迭代停止条件为||x1-x2||<w
double Power_method_1(vector<vector<double>> A, vector<double> &u);                                  //幂法迭代一次
double Power_method(vector<vector<double>> A, vector<double>& u,double w,int &k);                           //幂法求模最大特征值和特征向量
vector<vector<double>> up_Hessenberg_with_household(vector<vector<double>> A);                       //Household变换法求上Hessenberg分解（算法6.4.1）次对角线要非零
vector<vector<double>> QR_iteration(vector<vector<double>>H, vector<vector<double>> &P);//算法6.4.2双重步位移的QR迭代,返回变换后的H,保存P
void m_l(vector<vector<double>> A, int& m, int& l);                        //算法6.4.3确定最大的m和l
vector<vector<double>> Implicit_QR_algo(vector<vector<double>> A, double u,int &k);//隐式QR算法求解A的特征值,机器精度为u
void cout_diag_characteristic_value(vector<vector<double>> A);//输出拟对角阵的特征值





