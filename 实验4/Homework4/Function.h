#pragma once
#include<iostream>
#include <vector>
#include<math.h>
#define double long double
using namespace std;

vector<vector<double>> Matrix_add(vector<vector<double>> a, vector<vector<double>> b);                 //�������
vector<vector<double>> Matrix_Mul(vector<vector<double>> a, vector<vector<double>> b);                 //�������(Ĭ�Ͽ�����ˣ����������жϣ�
vector<double> Matrix_Mul_Vector(vector<vector<double>> A, vector<double> x);            //���������
vector<double> Vector_Mul_Matrix(vector<double> x ,vector<vector<double>> A);            //�����˾���
vector<vector<double>> LVector_Mul_HVector(vector<double> x, vector<double> y);          //��������������

vector<vector<double>> Matrix_T(vector<vector<double>> a);                  //�����ת��
vector<vector<double>> Matrix_Mul_num(vector<vector<double>> a, double b);               //��������
vector<vector<double>> Vector_Mul_num(vector<vector<double>> a, double b);               //��������
vector<vector<double>> Matrix_Divid_num(vector<vector<double>> a, double b);                           //�������һ����
vector<vector<double>> Matrix_Subtracition(vector<vector<double>> a, vector<vector<double>> b);   //�������,a-b
vector<double> Vector_Subtracition(vector<double> a, vector<double> b);                          //�������
vector<double> Vector_Add(vector<double> a, vector<double> b);                                  //�������
vector<vector<double>> make_zero_matrix(int a, int b);                                                 //����0����
vector<vector<double>> I(int b);                                                                       //����b�׵�λ����
double det(vector<vector<double>> A);                                                            //��Ϊ�����Ǿ�����det��A��
vector<vector<double>> read_Matrix(int m, int n);                                                      // ����m * n���󣨣�
vector<vector<double>> Matrix_block(vector<vector<double>> A, int start_x, int start_y, int block_x, int block_y);//�ֿ����(��start_x,start_y��ʼ����СΪblock_x * block_y�ķֿ����)
vector<double> Upper_triangle(vector<vector<double>> U, vector<double> b);                       //����������󷽳���Ux = b
vector<double> Lower_triangle(vector<vector<double>> L, vector<double> b);                       //����������󷽳���Lx = b
void cout_Matrix(vector<vector<double>> A);                                                //�������
void cout_Lvector(vector<double> b);                                                            //���������
void cout_Hvector(vector<double> b);                                                            //���������
void swap_row(vector<vector<double>>& A, int k, int p);                                  //����m * n�׷����k�к͵�p��
void swap_col(vector<vector<double>>& A, int k, int p);                                  //����m * n�׷����k�к͵�p��
void swap_vector(vector<double>& b, int k, int p);                                                     //��������b�ĵ�k,p��Ԫ��
void max_subMatrix(vector<vector<double>> A, int j, int t, int m, int n, int& p, int& q);              //��ĳ���Ӿ���ij--mn�������ֵ�����Ԫ��λ��
void getLU(vector<vector<double>> A, vector<vector<double>>& L, vector<vector<double>>& U);//ȡL,U����
vector<double> Gauss_no(vector<vector<double>> A,vector<double> b);                       //��ѡ��Ԫ��Gauss��ȥ��
vector<double> Gauss_all(vector<vector<double>> A,  vector<double> b);                    //ȫ��ԪGauss��ȥ��
vector<double> Gauss_col(vector<vector<double>> A, vector<double> b);                     //����ԪGauss��ȥ��
//void Gauss_col_LU(vector<vector<double>> A, const int& n, vector<double>& L, vector<double>& U);     //����Ԫ���Ƿֽ���L��U
vector<double> Square_LLT(vector<vector<double>> A, vector<double> b);                    //ƽ������ȥ��
vector<double> Square_LDLT(vector<vector<double>> A,  vector<double> b);                  //�Ľ���ƽ������ȥ��
double Vector_norm( vector<double> b,double k);                                                 //��������
double Matrix_norm(vector<vector<double>> A,int k);                                      //������
double Matrix_norm_1_pro(vector<vector<double>> A);                                              //�Ż������ƾ���A����������
vector<vector<double>> sign_Matrix(vector<vector<double>> w);                            //�жϾ���Ԫ�ط���
double K(vector<vector<double>> A,int k);                                                       //��A��������  
void House(vector<double> x, vector<double> &v, double &beta);                                  //����Householder�任
void givens(double c, double s, double a, double b);                                                   //����givens�任
vector<vector<double>> Q_R(vector<vector<double>> A, vector<double>& d);                 //����A��QR�ֽ�
vector<double> Q_R_solve(vector<vector<double>> A, vector<double> b);                    //��QR�ֽ�ⷽ��
vector<double> Jacobi_Iterative_method(vector<vector<double>> A, vector<double> b, vector<double> x);            //Jacobi�������ⷽ��
vector<double> G_S_Iterative_method(vector<vector<double>> A, vector<double> b, vector<double> x);               //G-S�������ⷽ��
vector<double> SOR_Iterative_method(vector<vector<double>> A, vector<double> b, vector<double> x, double w);     //SOR�������ⷽ��
vector<vector<double>> Jacobi_Iterative_method2(vector<vector<double>> U, vector<vector<double>> F, vector<vector<double>> G);          //Jacobi�������ⷽ�̣���άƫ΢�ַ���
vector<vector<double>> G_S_Iterative_method2(vector<vector<double>> A, vector<vector<double>> F, vector<vector<double>> G);               //G-S�������ⷽ��,��άƫ΢�ַ���
vector<vector<double>> SOR_Iterative_method2(vector<vector<double>> A, vector<vector<double>> F, vector<vector<double>> G, double w);   //SOR�������ⷽ��,��άƫ΢�ַ���
