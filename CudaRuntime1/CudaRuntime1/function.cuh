#include<stdio.h>
#include<stdlib.h>
#include<cuda_runtime.h>
#include<device_launch_parameters.h>
#include<math.h>
#include "cublas_v2.h"//动态链接库添加   cublas.lib;
#include <cusparse.h>//动态链接库添加    cusparse.lib; cublas.lib;
#include<time.h>

#define E 2.06e11					  //弹性模量
#define NODE 1701				      //结点数
#define ELE  1600					  //单元数
#define grid_length  2.5			 //网格尺寸
#define delta (3.03*grid_length)     //近场半径
#define s_crit 0.01                 //临界伸长率
#define L  15                        //长度尺度参数
#define threshold  1e-17                //CG迭代残差值
//以下几个数可先调大，运行算例，根据机器内存大小和算例网格情况选择
#define skip_size 50               //相关单元数组，也就是每个单元的相关单元个数最大值
#define skip_size_half 50			 //一半的相关单元数组
#define node_to_node_skip 80       //相关节点的数组大小，也就是每个结点的相关结点个数最大值

//此份代码的操作步骤：
//1.将网格文件拷贝至“mesh”文件夹中
//2.修改“function.cuh”文件中的参数：E、NODE、ELE、grid_length、delta等
//3.添加预制裂纹或屏蔽边界单元关系：“function.cu”中的“device_find_related_element”函数
//4.修改边界条件：“function.cu”中的“kernel_deal_bc”函数
//5.更新断键时屏蔽边界条件：“function.cu”中的“kernel_update_bond”函数
//6.输出反力F的范围修改：”function.cu“中的”export_result_load“反力范围修改
//7.“main.cu”中修改两个“while”循环中的”displacenment“值，添加加载步与加载条件
//8.运行即可

__host__ void host_read_FEM_information(double*node, int*element);

__host__ double host_cal_c0(double elatic, double delta_f);

__host__ void host_find_relevant_element(double*node, int*element, int*relevant_element,float* time);

__host__ void host_calculate_PeriFEM_stiffness(double*node, int* element, int*relevant_element, double* k_PeriFEM_stiffness,float* time);

__host__ void export_result_end(double*node, int*element, int* relevant_element, int*PeriFEM_bond,double*ans_CSR);

__host__ void host_cal_node_to_node_relation(int*element, int*relevant_element, int* node_to_node, int*CSR_value_size);

__host__ void host_cal_CSR_count_col(int CSR_value_size, int *k_all_CSR_col, int *k_all_CSR_count, int *node_to_node);

__host__ void update_bond(double*node, int*element, int*relevant_element, double*ans_CSR, int*PeriFEM_bond, int*PeriFEM_ele_broken, int *broken_ele_num,float*time);

__host__ void update_k_all(double*node, int*element, int*relevant_element, int*k_all_CSR_count,int* k_all_CSR_col,double* k_all_CSR_value, 
	int*PeriFEM_bond_backups,int* PeriFEM_bond,int*PeriFEM_ele_broken, int*broken_ele_num);

__host__ void host_bond_initialize(int num, int*mar);

__host__ void copy(double*k_all_CSR_value_ans, double*k_all_CSR_value,int CSR_value_size);

__host__ void Element_Pointer_Matrix(int*ele_1, int*ele_2, int*ele_P);//单元指示矩阵

__host__ void copy_bond(int*PeriFEM_bond_ans, int*PeriFEM_bond, int num);

__host__ void host_k_all(int*element, int*relevant_element, int*k_all_CSR_count, int*k_all_CSR_col, double*k_all_CSR_value, double * k_PeriFEM_stiffness);

__host__ void export_result_load(double*node, int*element, int* relevant_element, int*PeriFEM_bond, double*ans_CSR,double*F);

__host__ int Order_to_add(int num);

__host__ void host_cal_F(int*k_all_CSR_count, int*k_all_CSR_col, double*k_all_CSR_value, int CSR_value_size, double*ans_CSR, double*F);

__host__ void deal_bc(double*k_all_CSR_value, int*k_all_CSR_count, int*k_all_CSR_col, int CSR_value_size, double*R, double* node, double displacement,float* time_deal_bc);

__host__ void CG_iterative(double*k_all_CSR_value, int* k_all_CSR_count, int* k_all_CSR_col, int CSR_value_size, double* ans_CSR, double* R,int*num);

__host__ void update_Calculate_PeriFEM_stiffness(double*node_1, double*node_2, double c0,double* k_1,int*bond_backups, int*bond);

__device__ void device_find_related_element(int tid_num, double*node_centroid, double* node_d, int* element_d, int*relevant);

__host__ __device__ void shape_function(double r, double s, double*N_1);

__host__ __device__ void max_min(double* temp, double a, double b);

__device__ void device_update_bond(double*node_1, double* node_1_cal, double* node_2, double* node_2_cal, int* bond,int*pointer);

__host__ __device__ void Jacobi(double*node_1, double r, double s, double*J_det);

__host__ __device__ void real_location(double*N_1, double*node_1, double*real_location_1);

__host__ __device__ void D_martrix(double*real_location_2, double c0,double*D);

__host__ __device__ void B_martrix(double*N_1, double*N_2, double* B);

__device__ void device_Calculate_PeriFEM_stiffness(double*node_1, double*node_2,double c0, double* k_1);

__host__ __device__ void zeros_d(int num, double*mar);

__host__ __device__ void zeros(int num, int*mar);

__host__ __device__ void negative(int num, int*mar);

__host__ __device__ void Sort(int num_f, int*mart);

__global__ void kernel_find_relevant_element(double*node_d, int* element_d, int* relevant_element_d);

__global__ void kernel_Calculate_PeriFEM_stiffness(double*node_d, int* element_d, int* relevant_element_d,double c0, double* k_PeriFEM_stiffness_d);

__global__ void kernel_update_bond(double*node_d, int* element_d, int* relevant_element_d, double* ans_CSR_d, int* PeriFEM_bond_d, int*PeriFEM_ele_broken_d);

__global__ void kernel_deal_bc(double*k_all_CSR_value, int*k_all_CSR_count, int*k_all_CSR_col, double*R, double* node, double displacement);

__global__ void kernel_damage (int*element_d,double* damage_d, double*ele_bond_d);