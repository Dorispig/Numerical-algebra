#pragma once
#include<iostream>
#include <vector>
#include<math.h>
#define double long double
using namespace std;

vector<vector<double>> Matrix_add(vector<vector<double>> a, vector<vector<double>> b);                 //�������
vector<vector<double>> Matrix_Mul(vector<vector<double>> a, vector<vector<double>> b,const int& m, const int& n, const int& p); //�������(Ĭ�Ͽ�����ˣ����������жϣ�
vector<double> Matrix_Mul_Vector(vector<vector<double>> A, vector<double> x, int m, int n);            //���������
vector<double> Vector_Mul_Matrix(vector<double> x ,vector<vector<double>> A, int m, int n);            //�����˾���
vector<vector<double>> LVector_Mul_HVector(vector<double> x, vector<double> y, int m, int n);    //��������������

vector<vector<double>> Matrix_T(vector<vector<double>> a,const int& m, const int& n);                  //�����ת��
vector<vector<double>> Matrix_Mul_num(vector<vector<double>> a, double b);                             //��������
vector<vector<double>> Matrix_Divid_num(vector<vector<double>> a, double b);                           //�������һ����
vector<vector<double>> Matrix_Subtracition(vector<vector<double>> a, vector<vector<double>> b,int m,int n);        //�������,a-b
vector<double> Vector_Subtracition(vector<double> a, vector<double> b,int m);                          //�������
vector<vector<double>> make_zero_matrix(int a, int b);                                                 //����0����
vector<vector<double>> I(int b);                                                                         //����b�׵�λ����
double det(vector<vector<double>> A,int n);                                                            //��Ϊ�����Ǿ�����det��A��
vector<vector<double>> read_Matrix(int m, int n);                                                      // ����m * n���󣨣�
vector<vector<double>> Matrix_block(vector<vector<double>> A, int start_x, int start_y, int block_x, int block_y);//�ֿ����(��start_x,start_y��ʼ����СΪblock_x * block_y�ķֿ����)
vector<double> Upper_triangle(int n,vector<vector<double>> U, vector<double> b);                       //����������󷽳���Ux = b
vector<double> Lower_triangle(int n,vector<vector<double>> L, vector<double> b);                       //����������󷽳���Lx = b
void cout_Matrix(vector<vector<double>> A,int m,int n);                                                //�������
void cout_Lvector(vector<double> b, int n);                                                            //���������
void cout_Hvector(vector<double> b, int n);                                                            //���������
void swap_row(vector<vector<double>>& A, int m, int n, int k, int p);                                  //����m * n�׷����k�к͵�p��
void swap_col(vector<vector<double>>& A, int m, int n, int k, int p);                                  //����m * n�׷����k�к͵�p��
void swap_vector(vector<double>& b, int k, int p);                                                     //��������b�ĵ�k,p��Ԫ��
void max_subMatrix(vector<vector<double>> A, int j, int t, int m, int n, int& p, int& q);              //��ĳ���Ӿ���ij--mn�������ֵ�����Ԫ��λ��
void getLU(vector<vector<double>> A, const int& n, vector<vector<double>>& L, vector<vector<double>>& U);//ȡL,U����
vector<double> Gauss_no(vector<vector<double>> A,const int &n,vector<double> b);                       //��ѡ��Ԫ��Gauss��ȥ��
vector<double> Gauss_all(vector<vector<double>> A, const int& n, vector<double> b);                    //ȫ��ԪGauss��ȥ��
vector<double> Gauss_col(vector<vector<double>> A, const int& n,vector<double> b);                     //����ԪGauss��ȥ��
//void Gauss_col_LU(vector<vector<double>> A, const int& n, vector<double>& L, vector<double>& U);       //����Ԫ���Ƿֽ���L��U
vector<double> Square_LLT(vector<vector<double>> A,const int& n, vector<double> b);                    //ƽ������ȥ��
vector<double> Square_LDLT(vector<vector<double>> A, const int& n, vector<double> b);                  //�Ľ���ƽ������ȥ��
double Vector_norm( vector<double> b, int n,double k);                                                 //��������
double Matrix_norm(vector<vector<double>> A, int m, int n,int k);                                      //������
double Matrix_norm_1_pro(vector<vector<double>> A,int n);                                              //�Ż������ƾ���A����������
vector<vector<double>> sign_Matrix(vector<vector<double>> w, int m, int n);                            //�жϾ���Ԫ�ط���
double K(vector<vector<double>> A, int n,int k);                                                       //��A��������  
void House(vector<double> x, vector<double> &v, double &beta, int n);                                  //����Householder�任
void givens(double c, double s, double a, double b);                                                   //����givens�任
vector<vector<double>> Q_R(vector<vector<double>> A, vector<double>& d, int m, int n);                         //����A��QR�ֽ�
vector<double> Q_R_solve(vector<vector<double>> A, vector<double> b, int m, int n);                      //��QR�ֽ�ⷽ��




