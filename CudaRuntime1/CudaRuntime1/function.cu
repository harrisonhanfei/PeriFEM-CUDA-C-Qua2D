#include"function.cuh"
__host__ void host_read_FEM_information(double*node, int*element)//读取结点与单元信息
{
	FILE*fp_node = fopen(".\\mesh\\node_infor.txt", "r");//.\\mesh\\node_infor.txt
	if (fp_node == NULL)
	{
		printf("读取结点信息失败！请检查该目录下是否有node_infor.txt文件！\n");
		return;
	}
	for (int i = 0; i < NODE; i++)
	{

		for (int j = 0; j < 3; j++)
		{
			fscanf(fp_node, "%lf", &node[j + i * 3]);
		}
		fscanf(fp_node, "\n");
	}
	fclose(fp_node);

	FILE*fp_element = fopen(".\\mesh\\ele_infor.txt", "r");
	if (fp_element == NULL)
	{
		printf("读取单元信息失败！请检查该目录下是否有ele_infor.txt文件！\n");
		return;
	}
	for (int i = 0; i < ELE; i++)
	{

		for (int j = 0; j < 5; j++)
		{
			fscanf(fp_element, "%d", &element[j + i * 5]);
		}
		fscanf(fp_element, "\n");
	}
	fclose(fp_element);
	printf("有限元节点信息输入完毕。\n");
}
__host__ __device__ void zeros_d(int num, double*mar)//数组置零
{
	for (int i = 0; i < num; i++)
	{
		mar[i] = 0.0;
	}
}
__host__ __device__ void zeros(int num, int*mar)//数组置零
{
	for (int i = 0; i < num; i++)
	{
		mar[i] = 0;
	}
}
__host__ __device__ void negative(int num, int*mar)//数组置-1
{
	for (int i = 0; i < num; i++)
	{
		mar[i] = -1;
	}
}
__host__ __device__ void Sort(int num_f, int*mart)//数组排序，又小到大
{
	for (int n = 0; n < num_f; n++)
	{
		for (int j = 0; j < num_f - 1 - n; j++)
		{
			if (mart[j] > mart[j + 1])
			{
				int b = mart[j + 1];
				mart[j + 1] = mart[j];
				mart[j] = b;
			}
		}
	}
}
__device__ void device_find_related_element(int tid_num,double*node_centroid, double* node_d, int* element_d, int*relevant)
{
	relevant[0] = 0;
	for (int tid = 0; tid < ELE; tid++)
	{
		double	node_centroid_1[2];
		zeros_d(2, node_centroid_1);
		for (int i = 0; i < 4; i++)
		{
			node_centroid_1[0] += node_d[1 + 3 * (element_d[i + 1 + tid * 5] - 1)];
			node_centroid_1[1] += node_d[2 + 3 * (element_d[i + 1 + tid * 5] - 1)];
		}
		node_centroid_1[0] /= 4;
		node_centroid_1[1] /= 4;
		double a_1 = node_centroid_1[0] - node_centroid[0];
		double a_2 = node_centroid_1[1] - node_centroid[1];
		double dis_1_2 = sqrt(a_1*a_1 + a_2 * a_2);
		//if (fabs(a_2) > 1e-8)//增加的裂纹或屏边边界单元关系
		//{
		//	double x, k_1, k_2;
		//	x = (0.5 - node_centroid[1])*a_1 / a_2 + node_centroid[0];
		//	k_1 = (node_centroid[1]-0.5) / (node_centroid[0] - 10.0);
		//	k_2 = (node_centroid_1[1] -0.5) / (node_centroid_1[0] - 10.0);
		//	if ((x <=0.5) && ((k_1 < 0 && k_2>0) || (k_2 < 0 && k_1>0)))continue;
		//}//origin
		//new_crack
		double crack0[4] = { 0.0, 50.0, 50.0, 50.0 };//裂纹线段俩端点的x1，y1，x2，y2
		double my_x[2], my_y[2], my_a[2], my_b[2];
		max_min(my_x, node_centroid[0], node_centroid_1[0]); //x0, x1
		max_min(my_y, node_centroid[1], node_centroid_1[1]); //y0, y1
		max_min(my_a, crack0[0], crack0[2]);
		max_min(my_b, crack0[1], crack0[3]);
		
		if (my_x[1] > my_a[0] || my_x[0] < my_a[1] || my_y[1] > my_b[0] || my_y[0] < my_b[1])
		{
			//printf("%d\r", tid_num);
		}
		else
		{
			//printf("1");
			double xy0_ab0, xy0_ab1, ab0_xy0, ab0_xy1;
			xy0_ab0 = (crack0[0] - node_centroid[0]) * (node_centroid_1[1] - node_centroid[1]) - (node_centroid_1[0] - node_centroid[0]) * (crack0[1] - node_centroid[1]);
			xy0_ab1 = (crack0[2] - node_centroid[0]) * (node_centroid_1[1] - node_centroid[1]) - (node_centroid_1[0] - node_centroid[0]) * (crack0[3] - node_centroid[1]);
			ab0_xy0 = (node_centroid[0] - crack0[0]) * (crack0[3] - crack0[1]) - (crack0[2] - crack0[0]) * (node_centroid[1] - crack0[1]);
			ab0_xy1 = (node_centroid_1[0] - crack0[0]) * (crack0[3] - crack0[1]) - (crack0[2] - crack0[0]) * (node_centroid_1[1] - crack0[1]);
			if (xy0_ab0 * xy0_ab1 < 0 && ab0_xy0 * ab0_xy1 < 0)
				continue;
		}
		if (dis_1_2 <= delta)
		{
			relevant[0]++;
			if (tid_num == tid)relevant[1] = relevant[0];
			relevant[relevant[0]+1] = tid;
		}
	}
}

__host__ __device__ void max_min(double* temp, double a, double b)//比较大小
{
	if (a > b)
	{
		temp[0] = a;
		temp[1] = b;
	}
	else
	{
		temp[0] = b;
		temp[1] = a;
	}
}

__host__ void host_bond_initialize(int num, int*mar)
{
	for (int i = 0; i < num; i++)
	{
		mar[i] = 1;
		if (i % 17 == 0)mar[i] = 16;
	}
}
__global__ void kernel_find_relevant_element(double*node_d, int* element_d, int* relevant_element_d)
{
	int tid = threadIdx.x + threadIdx.y * blockDim.x + blockIdx.x * blockDim.x * blockDim.y + blockIdx.y * gridDim.x * blockDim.x * blockDim.y;
	if (tid < ELE)
	{
		double	node_centroid[2];
		zeros_d(2, node_centroid);
		for (int i = 0; i < 4; i++)
		{
			node_centroid[0] += node_d[1 + 3 * (element_d[i + 1 + tid * 5] - 1)];
			node_centroid[1] += node_d[2 + 3 * (element_d[i + 1 + tid * 5] - 1)];
		}
		node_centroid[0] /= 4;
		node_centroid[1] /= 4;
		int relevant[skip_size];
		device_find_related_element(tid,node_centroid, node_d, element_d, relevant);
		for (int m = 0; m < relevant[0] + 2; m++)
		{
			relevant_element_d[m + skip_size * tid] = relevant[m];
		}
	}
}
__host__ void host_find_relevant_element(double*node, int*element, int*relevant_element,float* time)
{
	double*node_d;            cudaMalloc((void**)&node_d, 3 * NODE * sizeof(double));
	int*element_d;            cudaMalloc((void**)&element_d, 5 * ELE * sizeof(int));
	int*relevant_element_d;   cudaMalloc((void**)&relevant_element_d, ELE * skip_size * sizeof(int));

	cudaMemcpy(node_d, node, 3 * NODE * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(element_d, element, 5 * ELE * sizeof(int), cudaMemcpyHostToDevice);

	dim3 block(32, 32);   dim3 grid(32, 32);////二维网格   二维线程块

	cudaEvent_t start_cal_k, stop;
	cudaEventCreate(&start_cal_k);
	cudaEventCreate(&stop);
	cudaEventRecord(start_cal_k, 0);

	kernel_find_relevant_element << <grid, block >> > (node_d, element_d, relevant_element_d);//并行寻找相关单元
	
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float cal_rele_ele_time;
	cudaEventElapsedTime(&cal_rele_ele_time, start_cal_k, stop);//单位ms
	time[0] = cal_rele_ele_time;
	printf("GPU寻找相关单元用时: %f ms\n", cal_rele_ele_time);
	cudaEventDestroy(start_cal_k);
	cudaEventDestroy(stop);

	cudaMemcpy(relevant_element, relevant_element_d, ELE * skip_size * sizeof(int), cudaMemcpyDeviceToHost);
	cudaFree(node_d);
	cudaFree(element_d);
	cudaFree(relevant_element_d);
}
__host__ __device__ void shape_function(double r, double s, double*N_1)
{
	zeros_d(16, N_1);
	//printf("%lf\n", r);
	N_1[0] = N_1[9] = 1.0 / 4.0 * (1.0 - r) * (1.0 - s);
	N_1[2] = N_1[11] = 1.0 / 4.0 * (1.0 + r) * (1.0 - s);
	N_1[4] = N_1[13] = 1.0 / 4.0 * (1.0 + r) * (1.0 + s);
	N_1[6] = N_1[15] = 1.0 / 4.0 * (1.0 - r) * (1.0 + s);
	//printf("%lf\n", N_1[0]);
}
__host__ __device__ void Jacobi(double*node_1, double r, double s, double*J_det)
{
	double dsdr[8], J[4];

	dsdr[0] = (s - 1) / 4;
	dsdr[1] = (1 - s) / 4;
	dsdr[2] = (1 + s) / 4;
	dsdr[3] = (-1) * (1 + s) / 4;
	dsdr[4] = (r - 1) / 4;
	dsdr[5] = (-1) * (1 + r) / 4;
	dsdr[6] = (1 + r) / 4;
	dsdr[7] = (1 - r) / 4;

	zeros_d(4, J);//雅可比矩阵赋初值0

	for (int r = 0; r < 2; r++)//3层循环实现矩阵乘法 。求： 雅可比矩阵
	{
		for (int c = 0; c < 2; c++)
		{
			for (int n = 0; n < 4; n++)
			{
				J[2 * r + c] += dsdr[4 * r + n] * node_1[2 * n + c];
			}
		}
	}
	J_det[0] = J[0] * J[3] - J[1] * J[2];//雅可比矩阵的值
}
__host__ __device__ void real_location(double*N_1, double*node_1, double*real_location_1)
{
	zeros_d(2, real_location_1);
	for (int i = 0; i < 4; i++)
	{
		real_location_1[0] += N_1[2 * i] * node_1[2 * i];
		real_location_1[1] += N_1[2 * i + 9] * node_1[2 * i + 1];
	}
}
__host__ __device__ void D_martrix(double*real_location_2,double c0, double*D)
{
	double relative_distance = sqrt(real_location_2[0] * real_location_2[0] + real_location_2[1] * real_location_2[1]);
	if (relative_distance < 1e-8)
	{
		zeros_d(4, D);
		return;
	}
	double c_0 = 0.5 * c0*exp(-relative_distance * L / delta);
	D[0] = c_0 * real_location_2[0] * real_location_2[0];
	D[1] = c_0 * real_location_2[0] * real_location_2[1];
	D[2] = c_0 * real_location_2[0] * real_location_2[1];
	D[3] = c_0 * real_location_2[1] * real_location_2[1];
}
__host__ __device__ void B_martrix(double*N_1, double*N_2, double* B)
{
	zeros_d(32, B);
	for (int i = 0; i < 4; i++)
	{
		B[2 * i] = B[2 * i + 17] = N_1[2 * i];
		B[8 + 2 * i] = B[2 * i + 25] = -N_2[2 * i];
	}
}
__device__ void device_Calculate_PeriFEM_stiffness(double*node_1, double*node_2, double c0,double* k_1)
{
	double Gauss_Point[2] = { -0.5773502691896, 0.5773502691896 };//二阶高斯积分点
	//double Weight[2] = { 1.0,1.0 };//二阶高斯积分的权重比
	for (int ii = 0; ii < 2; ii++)
	{
		for (int ij = 0; ij < 2; ij++)
		{
			double N_1[16], J_det_1[1], real_location_1[2];
			shape_function(Gauss_Point[ii], Gauss_Point[ij], N_1);
			Jacobi(node_1, Gauss_Point[ii], Gauss_Point[ij], J_det_1);
			real_location(N_1, node_1, real_location_1);
			for (int ji = 0; ji < 2; ji++)
			{
				for (int jj = 0; jj < 2; jj++)
				{
					double N_2[16], J_det_2[1], real_location_2[2], D[4], B[32], B_t[32];
					shape_function(Gauss_Point[ji], Gauss_Point[jj], N_2);
					Jacobi(node_2, Gauss_Point[ji], Gauss_Point[jj], J_det_2);
					real_location(N_2, node_2, real_location_2);
					real_location_2[0] = real_location_2[0] - real_location_1[0];
					real_location_2[1] = real_location_2[1] - real_location_1[1];
					//printf("%e\t%e\t\n", real_location_2[0], real_location_2[1]);
					D_martrix(real_location_2, c0,D);
					//printf("%e\t%e\t%e\t%e\n", D[0], D[1], D[2], D[3]);
					//printf("N1:%lf\t%lf\t%lf\t%lf\n", N_1[0], N_1[2], N_1[4], N_1[6]);
					//printf("N2:%lf\t%lf\t%lf\t%lf\n", N_2[0], N_2[2], N_2[4], N_2[6]);
					B_martrix(N_1, N_2, B);
					//printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", B[0], B[2], B[4], B[6], B[8], B[10], B[12], B[14]);
					for (int r = 0; r < 16; r++)
					{
						for (int n = 0; n < 2; n++)
						{
							B_t[2 * r + n] = B[16 * n + r];
						}
					}
					double Middle[32];
					zeros_d(32, Middle);
					/*           计算B_t（16*2）矩阵 * D（2*2）矩阵          */
					for (int r = 0; r < 16; r++)
					{
						for (int c = 0; c < 2; c++)
						{
							for (int n = 0; n < 2; n++)
							{
								Middle[2 * r + c] += B_t[2 * r + n] * D[c + 2 * n];
							}
						}
					}
					for (int r = 0; r < 16; r++)
					{
						for (int c = 0; c < 16; c++)
						{
							for (int n = 0; n < 2; n++)
							{
								k_1[16 * r + c] += J_det_1[0] * J_det_2[0] * Middle[2 * r + n] * B[c + 16 * n];
							}
						}
					}
				}
			}
		}
	}
}
__host__ void update_bond(double*node, int*element, int*relevant_element, double*ans_CSR, int*PeriFEM_bond,int*PeriFEM_ele_broken,
	int *broken_ele_num,float*time)
{
	double*node_d;            cudaMalloc((void**)&node_d, 3 * NODE * sizeof(double));
	int*element_d;            cudaMalloc((void**)&element_d, 5 * ELE * sizeof(int));
	int*relevant_element_d;   cudaMalloc((void**)&relevant_element_d, ELE * skip_size * sizeof(int));
	double*ans_CSR_d;           cudaMalloc((void**)&ans_CSR_d, 2 * NODE * sizeof(double));
	int*PeriFEM_bond_d;       cudaMalloc((void**)&PeriFEM_bond_d, ELE * skip_size * 17 * sizeof(int));
	int*PeriFEM_ele_broken_d;   cudaMalloc((void**)&PeriFEM_ele_broken_d, ELE * skip_size  * sizeof(int));

	int*PeriFEM_ele_broken_inner = (int*)malloc(ELE * skip_size * sizeof(int));

	cudaMemcpy(node_d, node, 3 * NODE * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(element_d, element, 5 * ELE * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(relevant_element_d, relevant_element, ELE * skip_size * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ans_CSR_d, ans_CSR, 2 * NODE * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(PeriFEM_bond_d, PeriFEM_bond, ELE * skip_size * 17 * sizeof(int), cudaMemcpyHostToDevice);

	dim3 block(32, 32);   dim3 grid(32, 32);////二维网格   二维线程块
	cudaEvent_t start_cal_k, stop;
	cudaEventCreate(&start_cal_k);
	cudaEventCreate(&stop);
	cudaEventRecord(start_cal_k, 0);
	kernel_update_bond << <grid, block >> > (node_d, element_d, relevant_element_d, ans_CSR_d, PeriFEM_bond_d,PeriFEM_ele_broken_d);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float cal_rele_ele_time;
	cudaEventElapsedTime(&cal_rele_ele_time, start_cal_k, stop);//单位ms
	time[0] = cal_rele_ele_time;

	cudaEventDestroy(start_cal_k);
	cudaEventDestroy(stop);
	cudaMemcpy(PeriFEM_bond, PeriFEM_bond_d, ELE * skip_size * 17 * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(PeriFEM_ele_broken_inner, PeriFEM_ele_broken_d, ELE * skip_size * sizeof(int), cudaMemcpyDeviceToHost);
	int broken_num = 0;//损伤PeriFEM单元数
	for (int i = 0; i < ELE; i++)
	{
		if (PeriFEM_ele_broken_inner[1 + i * skip_size] != 0)
		{
			for (int j = 0; j < PeriFEM_ele_broken_inner[1 + i * skip_size]+2; j++)
			{
				PeriFEM_ele_broken[j + skip_size * broken_num] = PeriFEM_ele_broken_inner[j + skip_size * i];

			}
			broken_num++;
		}
	}
	broken_ele_num[0] = broken_num;
	//printf("\n%d\n", broken_num);
	free(PeriFEM_ele_broken_inner);
	cudaFree(node_d);
	cudaFree(element_d);
	cudaFree(relevant_element_d);
	cudaFree(ans_CSR_d);
	cudaFree(PeriFEM_bond_d);
	cudaFree(PeriFEM_ele_broken_d);
}
__host__ void update_k_all(double*node, int*element, int*relevant_element_d, int*k_all_CSR_count, int* k_all_CSR_col, double* k_all_CSR_value,
	int*PeriFEM_bond_backups, int* PeriFEM_bond,int*PeriFEM_ele_broken, int*broken_ele_num)
{
	double c0 = host_cal_c0(E, delta);//根据材料参数计算指数型微模量系数的c0
	for (int m = 0; m < broken_ele_num[0]; m++)//1
	{
		printf("更新进度：%.3f%%\r", 100.0* m / broken_ele_num[0]);
		int tid = PeriFEM_ele_broken[skip_size * m];

		double node_1[8];//主单元的坐标
		for (int i = 0; i < 4; i++)
		{
			node_1[2 * i] = node[1 + 3 * (element[i + 1 + tid * 5] - 1)];
			node_1[2 * i + 1] = node[2 + 3 * (element[i + 1 + tid * 5] - 1)];
			//printf("%lf   %lf\n", node_1[2 * i], node_1[2 * i + 1]);
		}

		int relevant_element[skip_size];//相关单元

		for (int i = 0; i < relevant_element_d[tid*skip_size]; i++)
		{
			relevant_element[i] = relevant_element_d[2 + i + tid * skip_size];
			//printf("%d\t", relevant_element[i]);
		}
		//printf("\n");
		int relevant_element_broken[skip_size];//损伤单元
		for (int i = 0; i < PeriFEM_ele_broken[1+m*skip_size]; i++)//PeriFEM_ele_broken[skip_size * m+1]   损伤单元个数
		{
			relevant_element_broken[i] = PeriFEM_ele_broken[2 + i + m * skip_size];

		}
		
		int mark_broken_element[skip_size];//损伤标记
		for (int i = 0; i < PeriFEM_ele_broken[skip_size * m + 1]; i++)
		{
			int index;
			for (int j = 0; j < relevant_element_d[tid*skip_size]; j++)
			{
				if(relevant_element[j]== relevant_element_broken[i])
				{
					index = j;
					break;
				}
			}
			mark_broken_element[i] = index;
		}

		for (int i = 0; i < PeriFEM_ele_broken[skip_size * m + 1]; i++)//1   单元里坏了几个单元
		{
			double node_2[8];
			for (int j = 0; j < 4; j++)
			{
				node_2[2 * j] = node[1 + 3 * (element[j + 1 + relevant_element_broken[i] * 5] - 1)];
				node_2[2 * j + 1] = node[2 + 3 * (element[j + 1 + relevant_element_broken[i] * 5] - 1)];
			}

			int bond_backups[16];
			int bond[16];
			for (int j = 0; j < 16; j++)
			{
				bond_backups[j] = PeriFEM_bond_backups[1 + j + mark_broken_element[i] * 17 + tid * 17 * skip_size];
				bond[j] = PeriFEM_bond[1 + j + mark_broken_element[i] * 17 + tid * 17 * skip_size];
			}

			double k_1[16 * 16];
			zeros_d(16 * 16, k_1);
			update_Calculate_PeriFEM_stiffness(node_1, node_2,c0, k_1,bond_backups, bond);
			
			//更新总纲

			int ele_1[4], ele_2[4], node_point[16];//提出的结点都为真实结点编号（基于1的索引

			for (int j = 0; j < 4; j++)
			{
				ele_1[j] = element[j + 1 + 5 * tid];
				ele_2[j] = element[j + 1 + 5 * relevant_element_broken[i]];
			}

			Element_Pointer_Matrix(ele_1, ele_2, node_point);

			for (int ii = 0; ii < 16; ii++)
			{
				for (int j = 0; j < 16; j++)
				{
					int x = node_point[ii], y = node_point[j];
					int index;
					for (int n = k_all_CSR_count[x]; ; n++)
					{
						if (k_all_CSR_col[n] == y)
						{
							index = n;
							break;
						}
					}
					k_all_CSR_value[index] -= k_1[j + ii * 16];
				}
			}	
		}
	}
}
__host__ void update_Calculate_PeriFEM_stiffness(double*node_1, double*node_2,double c0, double* k_1, 
		 int*bond_backups, int*bond)
{
	double Gauss_Point[2] = { -0.5773502691896, 0.5773502691896 };//二阶高斯积分点
	for (int ii = 0; ii < 2; ii++)
	{
		for (int ij = 0; ij < 2; ij++)
		{
			double N_1[16], J_det_1[1], real_location_1[2];
			shape_function(Gauss_Point[ii], Gauss_Point[ij], N_1);
			Jacobi(node_1, Gauss_Point[ii], Gauss_Point[ij], J_det_1);
			real_location(N_1, node_1, real_location_1);
			for (int ji = 0; ji < 2; ji++)
			{
				for (int jj = 0; jj < 2; jj++)
				{
					if ( bond[jj + 2 * ji + 4 * ij + 8 * ii] ==bond_backups[jj + 2 * ji + 4 * ij + 8 * ii])continue;
					double N_2[16], J_det_2[1], real_location_2[2], D[4], B[32], B_t[32];
					shape_function(Gauss_Point[ji], Gauss_Point[jj], N_2);
					Jacobi(node_2, Gauss_Point[ji], Gauss_Point[jj], J_det_2);
					real_location(N_2, node_2, real_location_2);
					real_location_2[0] = real_location_2[0] - real_location_1[0];
					real_location_2[1] = real_location_2[1] - real_location_1[1];
					D_martrix(real_location_2, c0,D);
					B_martrix(N_1, N_2, B);
					for (int r = 0; r < 16; r++)
					{
						for (int n = 0; n < 2; n++)
						{
							B_t[2 * r + n] = B[16 * n + r];
						}
					}
					double Middle[32];
					zeros_d(32, Middle);
					for (int r = 0; r < 16; r++)
					{
						for (int c = 0; c < 2; c++)
						{
							for (int n = 0; n < 2; n++)
							{
								Middle[2 * r + c] += B_t[2 * r + n] * D[c + 2 * n];
							}
						}
					}
					for (int r = 0; r < 16; r++)
					{
						for (int c = 0; c < 16; c++)
						{
							for (int n = 0; n < 2; n++)
							{
								k_1[16 * r + c] += J_det_1[0] * J_det_2[0] * Middle[2 * r + n] * B[c + 16 * n];
							}
						}
					}
				}
			}
		}
	}
}



__host__ void copy(double*k_all_CSR_value_ans, double*k_all_CSR_value, int CSR_value_size)//复制数组
{
	for (int i = 0; i < CSR_value_size; i++)
	{
		k_all_CSR_value_ans[i] = k_all_CSR_value[i];
	}
}
__host__ void copy_bond(int*PeriFEM_bond_ans, int*PeriFEM_bond, int num)//复制数组
{
	for (int i = 0; i < num; i++)
	{
		PeriFEM_bond_ans[i] = PeriFEM_bond[i];
	}

}
__global__ void kernel_update_bond(double*node_d, int* element_d, int* relevant_element_d, double* ans_CSR_d, int* PeriFEM_bond_d, int*PeriFEM_ele_broken_d)
{
	int tid = threadIdx.x + threadIdx.y * blockDim.x + blockIdx.x * blockDim.x * blockDim.y + blockIdx.y * gridDim.x * blockDim.x * blockDim.y;
	if (tid < ELE)//屏蔽边界损伤处
	{
		double node_1[8];//主单元的原坐标
		for (int i = 0; i < 4; i++)
		{
			node_1[2 * i] = node_d[1 + 3 * (element_d[i + 1 + tid * 5] - 1)];
			node_1[2 * i + 1] = node_d[2 + 3 * (element_d[i + 1 + tid * 5] - 1)];
		}
		double node_1_cal[8];//主单元的计算后坐标
		for (int i = 0; i < 4; i++)
		{
			node_1_cal[2 * i] = ans_CSR_d[2 * (element_d[i + 1 + tid * 5] - 1)] + node_1[2 * i];
			node_1_cal[2 * i + 1] = ans_CSR_d[1 + 2 * (element_d[i + 1 + tid * 5] - 1)] + node_1[2 * i + 1];
		}
		int relevant_element[skip_size];
		for (int i = 0; i < relevant_element_d[tid*skip_size]; i++)
		{
			relevant_element[i] = relevant_element_d[2 + i + tid * skip_size];
		}
		int PeriFEM_ele_stat[skip_size];
		PeriFEM_ele_stat[0] = tid;
		PeriFEM_ele_stat[1] = 0;
		for (int i = 0; i < relevant_element_d[tid*skip_size]; i++)//1
		{
			double node_2[8];
			for (int j = 0; j < 4; j++)//提取从单元原坐标
			{
				node_2[2 * j] = node_d[1 + 3 * (element_d[j + 1 + relevant_element[i] * 5] - 1)];
				node_2[2 * j + 1] = node_d[2 + 3 * (element_d[j + 1 + relevant_element[i] * 5] - 1)];
			}
			double node_2_cal[8];
			for (int j = 0; j < 4; j++)
			{
				node_2_cal[2 * j] = ans_CSR_d[2 * (element_d[1 + j + relevant_element[i] * 5] - 1)] + node_2[2 * j];
				node_2_cal[2 * j + 1] = ans_CSR_d[1 + 2 * (element_d[1 + j + relevant_element[i] * 5] - 1)] + node_2[2 * j + 1];
			}
			int bond[17];
			for (int j = 0; j < 17; j++)
			{
				bond[j] = PeriFEM_bond_d[j + i * 17 + tid * 17 * skip_size];
			}
			int pointer[1] = { 0 };
			device_update_bond(node_1, node_1_cal, node_2, node_2_cal, bond, pointer);
			for (int j = 0; j < 17; j++)
			{
				PeriFEM_bond_d[j + i * 17 + tid * 17 * skip_size] = bond[j];
			}
			if (pointer[0]==1)
			{
				PeriFEM_ele_stat[1]++;
				PeriFEM_ele_stat[PeriFEM_ele_stat[1] + 1] = relevant_element[i];
			}
		}
		for (int i = 0; i < skip_size; i++)
		{
			PeriFEM_ele_broken_d[i + tid * skip_size] = PeriFEM_ele_stat[i];
		}
	}
}

__device__ void device_update_bond(double*node_1, double* node_1_cal, double* node_2, double* node_2_cal, int* bond, int*pointer)
{
	double Gauss_Point[2] = { -0.5773502691896, 0.5773502691896 };//二阶高斯积分点
	int broken_bond = 0;
	for (int ii = 0; ii < 2; ii++)
	{
		for (int ij = 0; ij < 2; ij++)
		{
			double N_1[16], real_location_1[2], real_location_1_cal[2];
			shape_function(Gauss_Point[ii], Gauss_Point[ij], N_1);
			real_location(N_1, node_1, real_location_1);
			real_location(N_1, node_1_cal, real_location_1_cal);
			for (int ji = 0; ji < 2; ji++)
			{
				for (int jj = 0; jj < 2; jj++)
				{
					double N_2[16], real_location_2[2], real_location_2_cal[2];
					shape_function(Gauss_Point[ji], Gauss_Point[jj], N_2);
					real_location(N_2, node_2, real_location_2);
					real_location(N_2, node_2_cal, real_location_2_cal);

					real_location_2[0] = real_location_2[0] - real_location_1[0];
					real_location_2[1] = real_location_2[1] - real_location_1[1];

					real_location_2_cal[0] = real_location_2_cal[0] - real_location_1_cal[0];
					real_location_2_cal[1] = real_location_2_cal[1] - real_location_1_cal[1];

					double dis_1_2 = sqrt(real_location_2[0] * real_location_2[0] + real_location_2[1] * real_location_2[1]);
					double dis_1_2_cal = sqrt(real_location_2_cal[0] * real_location_2_cal[0] + real_location_2_cal[1] * real_location_2_cal[1]);
					double s = (dis_1_2_cal - dis_1_2) / dis_1_2;
					if (s > s_crit && bond[1 + jj + ji * 2 + 4 * ij + 8 * ii] == 1)
					{
						pointer[0] = 1;
						broken_bond++;
						bond[1 + jj + ji * 2 + 4 * ij + 8 * ii] = 0;
					}
				}
			}
		}
	}
	bond[0] -= broken_bond;
}
__global__ void kernel_Calculate_PeriFEM_stiffness(double*node_d, int* element_d, int* relevant_element_d,double c0,double* k_PeriFEM_stiffness_d)
{
	int tid = threadIdx.x + threadIdx.y * blockDim.x + blockIdx.x * blockDim.x * blockDim.y + blockIdx.y * gridDim.x * blockDim.x * blockDim.y;
	if (tid < ELE)//
	{
		double node_1[8];//主单元的坐标
		for (int i = 0; i < 4; i++)
		{
			node_1[2 * i] = node_d[1 + 3 * (element_d[i + 1 + tid * 5] - 1)];
			node_1[2 * i + 1] = node_d[2 + 3 * (element_d[i + 1 + tid * 5] - 1)];
		}
		int relevant_element[skip_size_half];
		for (int i = 0; i < relevant_element_d[tid*skip_size+1]; i++)
		{
			relevant_element[i] = relevant_element_d[2 + i + tid * skip_size];
		}
		for (int i = 0; i < relevant_element_d[tid*skip_size+1]; i++)//1
		{
			double node_2[8];
			for (int j = 0; j < 4; j++)
			{
				node_2[2 * j] = node_d[1 + 3 * (element_d[j + 1 + relevant_element[i] * 5] - 1)];
				node_2[2 * j + 1] = node_d[2 + 3 * (element_d[j + 1 + relevant_element[i] * 5] - 1)];
			}
			double k_1[16 * 16];
			zeros_d(16 * 16, k_1);
			device_Calculate_PeriFEM_stiffness(node_1, node_2,c0, k_1);
			int count = 0;
			for (int m = 0; m < 16; m++)
			{
				for (int  n = 0; n < 16; n++)
				{
					k_PeriFEM_stiffness_d[count + 136 * i + tid * 136 * skip_size_half] = k_1[n+16*m];
					count++;
					if (m == n)break;
				}
			}
		}
	}
}
__host__ void host_calculate_PeriFEM_stiffness(double*node, int* element, int*relevant_element, double* k_PeriFEM_stiffness, float* time)
{
	double c0 = host_cal_c0(E, delta);//根据材料参数计算指数型微模量系数的c0

	double*node_d;                           cudaMalloc((void**)&node_d, 3 * NODE * sizeof(double));
	double*k_PeriFEM_stiffness_d;            cudaMalloc((void**)&k_PeriFEM_stiffness_d, ELE*skip_size_half * 136 * sizeof(double));
	int*element_d;            cudaMalloc((void**)&element_d, 5 * ELE * sizeof(int));
	int*relevant_element_d;   cudaMalloc((void**)&relevant_element_d, ELE * skip_size * sizeof(int));

	cudaMemcpy(node_d, node, 3 * NODE * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(element_d, element, 5 * ELE * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(relevant_element_d, relevant_element, ELE * skip_size * sizeof(int), cudaMemcpyHostToDevice);

	dim3 block(32, 32);   dim3 grid(32, 32);////二维网格   二维线程块
	cudaEvent_t start_cal_k, stop;
	cudaEventCreate(&start_cal_k);
	cudaEventCreate(&stop);
	cudaEventRecord(start_cal_k, 0);
	kernel_Calculate_PeriFEM_stiffness << <grid, block >> > (node_d, element_d, relevant_element_d,c0, k_PeriFEM_stiffness_d);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float cal_ele_k_time;
	cudaEventElapsedTime(&cal_ele_k_time, start_cal_k, stop);//单位ms
	time[0] = cal_ele_k_time;
	printf("GPU计算单刚用时: %f ms\n",  cal_ele_k_time);
	cudaEventDestroy(start_cal_k);
	cudaEventDestroy(stop);
	cudaMemcpy(k_PeriFEM_stiffness, k_PeriFEM_stiffness_d, ELE*skip_size_half * 136 * sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(node_d);
	cudaFree(element_d);
	cudaFree(k_PeriFEM_stiffness_d);
	cudaFree(relevant_element_d);
}
__host__ void host_cal_node_to_node_relation(int*element, int*relevant_element, int* node_to_node, int*CSR_value_size)
{
	int a = 0;
	zeros(NODE * node_to_node_skip, node_to_node);
	int max_node_occupy_element = 20;//一个结点占据的最大单元数,出错时调大此数
	int *node_element_show = (int*)malloc(sizeof(int) * max_node_occupy_element * NODE);
	int *count = (int*)malloc(sizeof(int)*NODE);
	zeros(max_node_occupy_element * NODE, node_element_show);
	zeros(NODE, count);//初始化

	for (int i = 1; i < ELE * 5; i += 5)
	{
		for (int j = 0; j < 4; j++)
		{
			node_element_show[(element[i + j] - 1) * max_node_occupy_element + count[element[i + j] - 1]] = (i / 5 + 1);
			count[element[i + j] - 1]++;
		}
	}
	
	int big_node_to_node = 0;
	for (int i = 0; i < NODE; i++)
	{
		int*extract_node = (int*)malloc(node_to_node_skip * sizeof(int));
		int*extract_element = (int*)malloc(max_node_occupy_element*skip_size * sizeof(int));
		negative(max_node_occupy_element*skip_size, extract_element);
		negative(node_to_node_skip, extract_node);
		
		int count_num = 0;
		for (int j = 0; j < count[i]; j++)
		{
			for (int m = 0; m < relevant_element[(node_element_show[j + i * max_node_occupy_element] - 1)*skip_size]; m++)
			{
				int sign = 0;
				for (int n = 0; n < 4 * skip_size; n++)
				{
					if (extract_element[n] == relevant_element[(node_element_show[j + i * max_node_occupy_element] - 1)*skip_size + m + 2])
					{
						sign = 1;
						break;
					}
				}
				if (sign == 0)
				{
					extract_element[count_num] = relevant_element[(node_element_show[j + i * max_node_occupy_element] - 1)*skip_size + m + 2];
					count_num++;
				}
			}
		}

		int count_num_1 = 0;
		for (int j = 0; j < count_num; j++)
		{
			for (int m = 0; m < 4; m++)
			{
				int sign = 0;
				for (int n = 0; n < node_to_node_skip; n++)
				{
					if (extract_node[n] == element[5 * (extract_element[j]) + m + 1])
					{
						sign = 1;
						break;
					}
				}
				if (sign == 0)
				{
					extract_node[count_num_1] = element[5 * (extract_element[j]) + m + 1];
					count_num_1++;
				}
			}
		}
		Sort(count_num_1, extract_node);
		a += count_num_1;
		for (int m = 0; m < count_num_1; m++)
		{
			node_to_node[m + node_to_node_skip * i] = extract_node[m];
		}
		if (count_num_1 > big_node_to_node)big_node_to_node = count_num_1;
		free(extract_node);
		free(extract_element);
	}
	CSR_value_size[0] = 4 * a;
	free(node_element_show);
	free(count);
	printf("big_node_to_node num: %d\n", big_node_to_node);
}
__host__ void host_cal_CSR_count_col(int CSR_value_size, int *k_all_CSR_col, int *k_all_CSR_count, int *node_to_node)
{
	k_all_CSR_count[0] = 0;
	zeros(CSR_value_size, k_all_CSR_col);

	for (int i = 0; i < NODE; i++)
	{
		int count_local;
		for (int n = 0; n < node_to_node_skip; n++)
		{
			if (node_to_node[n + node_to_node_skip * i] == 0)
			{
				count_local = n;
				break;
			}
		}
		k_all_CSR_count[2 * i + 1] = k_all_CSR_count[2 * i] + 2 * count_local;
		k_all_CSR_count[2 * i + 2] = k_all_CSR_count[2 * i + 1] + 2 * count_local;
	}//形成k_all_CSR_count行索引

	int a = 0;
	for (int i = 0; i < NODE; i++)//形成k_all_CSR_col列索引
	{
		int*extract = (int*)malloc(sizeof(int) * node_to_node_skip * 2);
		int count_local;
		for (int n = 0; n < node_to_node_skip; n++)
		{
			if (node_to_node[n + node_to_node_skip * i] == 0)
			{
				count_local = n;
				break;
			}
		}
		for (int n = 0; n < count_local; n++)
		{
			extract[2 * n] = 2 * node_to_node[n + node_to_node_skip * i] - 2;
			extract[2 * n + 1] = 2 * node_to_node[n + node_to_node_skip * i] - 1;
		}
		for (int m = 0; m < count_local * 2; m++)
		{
			k_all_CSR_col[m + a] = extract[m];
			k_all_CSR_col[m + count_local * 2 + a] = extract[m];
		}
		a += 4 * count_local;//步长
		free(extract);
	}
}
__host__ void Element_Pointer_Matrix(int*ele_1, int*ele_2, int*ele_P)//单元指示矩阵
{
	for (int i = 0; i < 4; i++)
	{
		ele_P[2 * i] = ele_1[i] * 2 - 2;
		ele_P[2 * i + 1] = ele_1[i] * 2 - 1;
	}
	for (int i = 4; i < 8; i++)
	{
		ele_P[2 * i] = ele_2[i - 4] * 2 - 2;
		ele_P[2 * i + 1] = ele_2[i - 4] * 2 - 1;
	}
}
__host__ int Order_to_add(int num)
{
	if (num==1)
	{
		return 0;
	}
	else
	{
		return num + Order_to_add(num - 1);
	}
}

__host__ void host_k_all(int*element, int*relevant_element, int*k_all_CSR_count, int*k_all_CSR_col, double*k_all_CSR_value, double * k_PeriFEM_stiffness)
{
	for (int m = 0; m < ELE; m++)//以单元做循环
	{
		printf("组装进度：%.3f%%\r", 100.0* m / (ELE));
		int* ele_1, *ele_2, *node_point;
		node_point = (int*)malloc(sizeof(int) * 16);
		ele_1 = (int*)malloc(sizeof(int) * 4);
		ele_2 = (int*)malloc(sizeof(int) * 4);//提出的结点都为真实结点编号（基于1的索引
		for (int i = 0; i < 4; i++)
		{
			ele_1[i] = element[i + 1 + 5 * m];
		}
		for (int i = 0; i < relevant_element[m*skip_size+1]; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				ele_2[j] = element[j + 1 + 5 * relevant_element[m*skip_size + 2 + i]];
			}
			Element_Pointer_Matrix(ele_1, ele_2, node_point);
			if (m == relevant_element[m*skip_size + 2 + i])
			{
				double *k_current = (double*)malloc(256 * sizeof(double));

				for (int ii = 0; ii < 16; ii++)
				{
					for (int jj = 0; jj < 16; jj++)
					{
						if (jj < ii)
						{
							k_current[jj + ii * 16] = k_PeriFEM_stiffness[Order_to_add(ii) + jj + 1 + i * 136 + m * skip_size_half * 136];
							k_current[ii + jj * 16] = k_PeriFEM_stiffness[Order_to_add(ii) + jj + 1 + i * 136 + m * skip_size_half * 136];
						}
						if (jj == ii)
						{
							k_current[ii + jj * 16] = k_PeriFEM_stiffness[Order_to_add(jj + 1) + i * 136 + m * skip_size_half * 136];
						}
					}
				}
				for (int ii = 0; ii < 16; ii++)
				{
					for (int j = 0; j < 16; j++)
					{
						int x = node_point[ii], y = node_point[j];
						int index;
						for (int n = k_all_CSR_count[x]; ; n++)
						{
							if (k_all_CSR_col[n] == y)
							{
								index = n;
								break;
							}
						}
						k_all_CSR_value[index] += k_current[j + ii * 16];
					}
				}
				free(k_current);
			}
			else
			{
				double *k_current = (double*)malloc(256 * sizeof(double));

				for (int ii = 0; ii < 16; ii++)
				{
					for (int jj = 0; jj < 16; jj++)
					{
						if (jj < ii)
						{
							k_current[jj + ii * 16] = 2.0*k_PeriFEM_stiffness[Order_to_add(ii) + jj + 1 + i * 136 + m * skip_size_half * 136];
							k_current[ii + jj * 16] = 2.0*k_PeriFEM_stiffness[Order_to_add(ii) + jj + 1 + i * 136 + m * skip_size_half * 136];
						}
						if (jj == ii)
						{
							k_current[ii + jj * 16] = 2.0*k_PeriFEM_stiffness[Order_to_add(jj + 1) + i * 136 + m * skip_size_half * 136];
						}
					}
				}
				for (int ii = 0; ii < 16; ii++)
				{
					for (int j = 0; j < 16; j++)
					{
						int x = node_point[ii], y = node_point[j];
						int index;
						for (int n = k_all_CSR_count[x]; ; n++)
						{
							if (k_all_CSR_col[n] == y)
							{
								index = n;
								break;
							}
						}
						k_all_CSR_value[index] += k_current[j + ii * 16];
					}
				}
				free(k_current);
			}

		}
		free(ele_1);
		free(ele_2);
		free(node_point);
	}
}
__host__ void deal_bc(double*k_all_CSR_value, int*k_all_CSR_count, int*k_all_CSR_col, int CSR_value_size, double*R, double* node, double displacement,float *time)
{
	double *value_d, *node_d, *R_d;
	int*count_d, *col_d;
	cudaMalloc((void**)&value_d, sizeof(double) *CSR_value_size);
	cudaMalloc((void**)&count_d, sizeof(int) * (NODE * 2 + 1));
	cudaMalloc((void**)&col_d, sizeof(int) * CSR_value_size);
	cudaMalloc((void**)&node_d, sizeof(double) * NODE * 3);
	cudaMalloc((void**)&R_d, sizeof(double)*NODE * 2);

	cudaMemcpy(value_d, k_all_CSR_value, sizeof(double) * CSR_value_size, cudaMemcpyHostToDevice);
	cudaMemcpy(count_d, k_all_CSR_count, sizeof(int) * (NODE * 2 + 1), cudaMemcpyHostToDevice);
	cudaMemcpy(col_d, k_all_CSR_col, sizeof(int) * CSR_value_size, cudaMemcpyHostToDevice);
	cudaMemcpy(node_d, node, sizeof(double) * NODE * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(R_d, R, sizeof(double) * NODE * 2, cudaMemcpyHostToDevice);
	cudaEvent_t start_cal_k, stop;
	cudaEventCreate(&start_cal_k);
	cudaEventCreate(&stop);
	cudaEventRecord(start_cal_k, 0);
	dim3 block(32, 32);//二维线程块
	dim3 grid(32, 32);//二维网格
	kernel_deal_bc << <grid, block >> > (value_d, count_d, col_d, R_d, node_d, displacement);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float cal_rele_ele_time;
	cudaEventElapsedTime(&cal_rele_ele_time, start_cal_k, stop);//单位ms
	time[0] = cal_rele_ele_time;

	cudaEventDestroy(start_cal_k);
	cudaEventDestroy(stop);

	cudaMemcpy(R, R_d, sizeof(double)*NODE * 2, cudaMemcpyDeviceToHost);
	cudaMemcpy(k_all_CSR_value, value_d, sizeof(double) * CSR_value_size, cudaMemcpyDeviceToHost);
	cudaFree(value_d);
	cudaFree(count_d);
	cudaFree(col_d);
	cudaFree(node_d);
	cudaFree(R_d);

}
__global__ void kernel_deal_bc(double*k_all_CSR_value, int*k_all_CSR_count, int*k_all_CSR_col, double*R, double* node, double displacement)
{
	int i = threadIdx.x + threadIdx.y * blockDim.x + blockIdx.x * blockDim.x * blockDim.y + blockIdx.y * gridDim.x * blockDim.x * blockDim.y;
	int index = 0;
	if (i < NODE)//
	{
		if (node[2 + i * 3] <  2.5 * grid_length)
		{	
			for (int m = k_all_CSR_count[2 * i+1]; ; m++)
			{
				if (k_all_CSR_col[m] == 2 * i+1)
				{
					index = m;
					break;
				}
			}
			k_all_CSR_value[index] *= 1e15;
			R[2 * i + 1] = - k_all_CSR_value[index] * displacement;
		}

		if (node[2 + i * 3] > 100 - 2.5 * grid_length)
		{
			for (int m = k_all_CSR_count[2 * i + 1]; ; m++)
			{
				if (k_all_CSR_col[m] == 2 * i+1)
				{
					index = m;
					break;
				}
			}
			k_all_CSR_value[index] *= 1e15;
			R[2 * i + 1] = k_all_CSR_value[index] * displacement;
		}

	}
}
__host__ void CG_iterative(double*k_all_CSR_value, int* k_all_CSR_count, int* k_all_CSR_col, int CSR_value_size, double* ans_CSR, double* R,int*num)
{
	double alpha = 1.0;
	double beat = 0.0;

	int size = NODE * 2;

	double *value_d, *x_d, *y_d, *ans_CSR_d;
	cudaMalloc((void**)&value_d, sizeof(double) *CSR_value_size);
	cudaMalloc((void**)&x_d, sizeof(double) * size);
	cudaMalloc((void**)&y_d, sizeof(double) * size);
	cudaMalloc((void**)&ans_CSR_d, sizeof(double) * size);

	cudaMemcpy(value_d, k_all_CSR_value, sizeof(double) * CSR_value_size, cudaMemcpyHostToDevice);
	cudaMemcpy(x_d, R, sizeof(double) * size, cudaMemcpyHostToDevice);

	cudaMemcpy(ans_CSR_d, ans_CSR, sizeof(double) * size, cudaMemcpyHostToDevice);
	int*count_d, *col_d;
	cudaMalloc((void**)&count_d, sizeof(int) * (size + 1));
	cudaMalloc((void**)&col_d, sizeof(int) * CSR_value_size);
	cudaMemcpy(count_d, k_all_CSR_count, sizeof(int) * (size + 1), cudaMemcpyHostToDevice);
	cudaMemcpy(col_d, k_all_CSR_col, sizeof(int) * CSR_value_size, cudaMemcpyHostToDevice);

	cusparseHandle_t     handle1 = NULL;
	cusparseSpMatDescr_t matA;
	cusparseDnVecDescr_t vecX, vecY;
	void*                dBuffer = NULL;
	size_t               bufferSize = 0;
	cusparseCreate(&handle1);// Create sparse matrix A in CSR format
	cusparseCreateCsr(&matA, size, size, CSR_value_size, count_d, col_d, value_d, CUSPARSE_INDEX_32I,
		CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);//创建cuspares里A的csr矩阵
	cusparseCreateDnVec(&vecX, size, x_d, CUDA_R_64F);
	// Create dense vector y
	cusparseCreateDnVec(&vecY, size, y_d, CUDA_R_64F);
	// allocate an external buffer if needed
	cusparseSpMV_bufferSize(handle1, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beat, vecY,
		CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize);
	//printf("%zd", bufferSize);//buffersize=41012    单位应该是 kb
	cudaMalloc(&dBuffer, bufferSize);

	cublasHandle_t handle;//cublas事件
	cublasCreate(&handle);//创建cublas事件

	double*r_0_d, *p_0_d, *r_1_d;
	cudaMalloc((void**)&r_0_d, sizeof(double) * size);
	cudaMalloc((void**)&p_0_d, sizeof(double) * size);
	cudaMalloc((void**)&r_1_d, sizeof(double) * size);

	cublasDcopy(handle, size, x_d, 1, r_0_d, 1);//复制x_d 到r_0_d
	cublasDcopy(handle, size, x_d, 1, p_0_d, 1);
	cublasDcopy(handle, size, x_d, 1, r_1_d, 1);

	double a_k, b_k, a_k_nagative;
	double	middle_a_0[1];
	double	middle_a_1[1];
	double	middle_b_0[1];
	double	middle_b_1[1];

	double  break_num[1];
	for (int ii = 0; ii < size; ii++)
	{
		printf("==迭代次数：%d==\r", ii);
		cublasDasum(handle, size, r_1_d, 1, break_num);//求向量内各元素绝对值和
		if (break_num[0] < threshold) { num[0] = ii; break; }//迭代停止
		cublasDdot(handle, size, r_0_d, 1, r_0_d, 1, middle_a_0);
		cusparseCreateDnVec(&vecX, size, p_0_d, CUDA_R_64F);
		cusparseSpMV(handle1, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beat, vecY,
			CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer);// vecY = alpha * matA*vecX + beat(0) * vecY
		cublasDdot(handle, size, p_0_d, 1, y_d, 1, middle_a_1);
		a_k = middle_a_0[0] / middle_a_1[0];
		a_k_nagative = -a_k;
		cublasDaxpy(handle, size, &a_k, p_0_d, 1, ans_CSR_d, 1);//ans_CSR_d=ans_CSR_d+alpha*p_0_d
		cublasDcopy(handle, size, r_0_d, 1, r_1_d, 1);
		cublasDaxpy(handle, size, &a_k_nagative, y_d, 1, r_0_d, 1);//r_0_d=r_0_d+a_k_nagative*middle_d
		cublasDdot(handle, size, r_1_d, 1, r_1_d, 1, middle_b_1);
		cublasDdot(handle, size, r_0_d, 1, r_0_d, 1, middle_b_0);
		b_k = middle_b_0[0] / middle_b_1[0];
		cublasDscal(handle, size, &b_k, p_0_d, 1);
		cublasDaxpy(handle, size, &alpha, r_0_d, 1, p_0_d, 1);
	}
	cublasDestroy(handle);
	cusparseDestroy(handle1);
	cudaMemcpy(ans_CSR, ans_CSR_d, sizeof(double) * size, cudaMemcpyDeviceToHost);//从设备端取回结果
	cudaFree(value_d);
	cudaFree(x_d);
	cudaFree(y_d);
	cudaFree(ans_CSR_d);
	cudaFree(count_d);
	cudaFree(col_d);
	cudaFree(r_0_d);
	cudaFree(p_0_d);
	cudaFree(r_1_d);
}
__host__ void export_result_end(double*node, int*element, int* relevant_element, int*PeriFEM_bond,double*ans_CSR)
{

	double*ele_bond = (double*)malloc(ELE * sizeof(double));
	for (int i = 0; i < ELE; i++)
	{
		int now_bond = 0;
		for (int j = 0; j < relevant_element[i*skip_size]; j++)
		{
			now_bond += PeriFEM_bond[17 * j + i * 17 * skip_size];
		}
		ele_bond[i] = double(now_bond) / 16.0 / relevant_element[i*skip_size];

	}

	double*damage = (double*)malloc(NODE * sizeof(double));
	
	int*element_d;            cudaMalloc((void**)&element_d, 5 * ELE * sizeof(int));
	double*damage_d;          cudaMalloc((void**)&damage_d, NODE * sizeof(double));
	double*ele_bond_d;        cudaMalloc((void**)&ele_bond_d, ELE * sizeof(double));

	cudaMemcpy(element_d, element, 5 * ELE * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ele_bond_d, ele_bond, ELE * sizeof(double), cudaMemcpyHostToDevice);

	dim3 block(32, 32);//二维线程块
	dim3 grid(32, 32);//二维网格
	
	kernel_damage << <grid, block >> > (element_d, damage_d, ele_bond_d);

	cudaMemcpy(damage, damage_d, NODE * sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(element_d);
	cudaFree(damage_d);
	cudaFree(ele_bond_d);
	free(ele_bond);

	FILE*fp;
	fp = fopen(".\\result\\GPU_date_end.dat", "w");
	fprintf(fp, "TITLE = Displacement_Field_Contour\nVARIABLES = X,Y,damage,Ux,Uy\nZONE N= %d E= %d ,F=FEPOINT, ET=QUADRILATERAL\n", NODE, ELE);
	for (int i = 0; i < NODE; i++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%e\t%e\n", node[3 * i + 1] + ans_CSR[2 * i], node[3 * i + 2] + ans_CSR[2 * i + 1], damage[i], ans_CSR[2 * i], ans_CSR[2 * i + 1]);
	}
	for (int i = 0; i < ELE; i++)
	{
		fprintf(fp, "%d %d %d %d\n", element[5 * i + 1], element[5 * i + 2], element[5 * i + 3], element[5 * i + 4]);
	}
	fclose(fp);
	free(damage);
}

__global__ void kernel_damage(int*element, double* damage, double*ele_bond)
{
	int tid = threadIdx.x + threadIdx.y * blockDim.x + blockIdx.x * blockDim.x * blockDim.y + blockIdx.y * gridDim.x * blockDim.x * blockDim.y;
	if (tid < NODE)
	{
		damage[tid] = 0.0;
		int ele_end[5];

		int count = 0;
		for (int e = 0; e < ELE; e++)
		{
			if ((tid + 1) == element[5 * e + 1] || (tid + 1) == element[5 * e + 2] || (tid + 1) == element[5 * e + 3] || (tid + 1) == element[5 * e + 4])
			{
				count++;
				ele_end[count - 1] = e;
			}
		}
		for (int n = 0; n < count; n++)
		{
			damage[tid] += ele_bond[ele_end[n]];
		}
		damage[tid] = 1.0 - damage[tid] / count;
	}

}
__host__ void export_result_load(double*node, int*element, int* relevant_element, int*PeriFEM_bond, double*ans_CSR,double*F)
{
	double*ele_bond = (double*)malloc(ELE * sizeof(double));
	for (int i = 0; i < ELE; i++)
	{
		int now_bond = 0;
		for (int j = 0; j < relevant_element[i*skip_size]; j++)
		{
			now_bond += PeriFEM_bond[17 * j + i * 17 * skip_size];
		}
		ele_bond[i] = double(now_bond) / 16.0 / relevant_element[i*skip_size];

	}

	double*damage = (double*)malloc(NODE * sizeof(double));

	int*element_d;            cudaMalloc((void**)&element_d, 5 * ELE * sizeof(int));
	double*damage_d;          cudaMalloc((void**)&damage_d, NODE * sizeof(double));
	double*ele_bond_d;        cudaMalloc((void**)&ele_bond_d, ELE * sizeof(double));

	cudaMemcpy(element_d, element, 5 * ELE * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ele_bond_d, ele_bond, ELE * sizeof(double), cudaMemcpyHostToDevice);

	dim3 block(32, 32);//二维线程块
	dim3 grid(32, 32);//二维网格

	kernel_damage << <grid, block >> > (element_d, damage_d, ele_bond_d);

	cudaMemcpy(damage, damage_d, NODE * sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(element_d);
	cudaFree(damage_d);
	cudaFree(ele_bond_d);
	free(ele_bond);

	FILE*fp;
	fp = fopen(".\\result\\GPU_damage_load.dat", "a");
	fprintf(fp, "TITLE = Displacement_Field_Contour\nVARIABLES = X,Y,damage,Ux,Uy,Fx,Fy\nZONE N= %d E= %d ,F=FEPOINT, ET=QUADRILATERAL\n", NODE, ELE);
	for (int i = 0; i < NODE; i++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", node[3 * i + 1], node[3 * i + 2], damage[i], ans_CSR[2 * i], ans_CSR[2 * i + 1], F[2 * i], F[2 * i + 1]);
	}
	for (int i = 0; i < ELE; i++)
	{
		fprintf(fp, "%d %d %d %d\n", element[5 * i + 1], element[5 * i + 2], element[5 * i + 3], element[5 * i + 4]);
	}
	fclose(fp);
	free(damage);

	//输出支反力
	double Fx = 0;
	double Fy = 0;
	int index;
	for (int i = 0; i < NODE; i++)
	{
		if ( node[2 + i * 3] > 100 - 2 * grid_length)
		{
			Fx += F[2 * i];
			Fy += F[2 * i + 1];
			index = i;
		}
	}
	FILE*fp_F;
	fp_F = fopen(".\\result\\GPU_F_load.dat", "a");
	fprintf(fp_F, "%lf\t%lf\t%lf\n", fabs(ans_CSR[2 * index + 1]), fabs(Fx), fabs(Fy));
	fclose(fp_F);
}

__host__ double host_cal_c0(double elatic, double delta_f)
{
	double c0 = 0;
	double h = delta_f / 100000;
	//printf("%lf,%lf,%lf\n", elatic, delta_f, h);
	for (int i = 0; i < 100000; i++)
	{
		c0 += h * (i*h + h / 2)* (i*h + h / 2)* (i*h + h / 2)* (i*h + h / 2)* (i*h + h / 2)*exp(-L * (i*h + h / 2) / delta_f);
	}
	c0 = 3.0*E / 3.1415926 / c0;
	return c0;
}

__host__ void host_cal_F(int*k_all_CSR_count, int*k_all_CSR_col, double*k_all_CSR_value,int CSR_value_size, double*ans_CSR, double*F)
{
	int size = NODE * 2;
	double alpha = 1.0;
	double beta = 0.0;
	
	// Device memory management
	double *value_d, *ans_CSR_d,*F_d;
	int*count_d, *col_d;
	cudaMalloc((void**)&value_d, sizeof(double) *CSR_value_size);
	cudaMalloc((void**)&ans_CSR_d, sizeof(double) * size);
	cudaMalloc((void**)&F_d, sizeof(double) * size);
	cudaMalloc((void**)&count_d, sizeof(int) * (size + 1));
	cudaMalloc((void**)&col_d, sizeof(int) * CSR_value_size);

	cudaMemcpy(value_d, k_all_CSR_value, sizeof(double) * CSR_value_size, cudaMemcpyHostToDevice);
	cudaMemcpy(ans_CSR_d, ans_CSR, sizeof(double) * size, cudaMemcpyHostToDevice);
	cudaMemcpy(count_d, k_all_CSR_count, sizeof(int) * (size + 1), cudaMemcpyHostToDevice);
	cudaMemcpy(col_d, k_all_CSR_col, sizeof(int) * CSR_value_size, cudaMemcpyHostToDevice);

	cusparseHandle_t     handle = NULL;
	cusparseSpMatDescr_t matA;
	cusparseDnVecDescr_t vecX, vecY;
	void*                dBuffer = NULL;
	size_t               bufferSize = 0;
	cusparseCreate(&handle);
		// Create sparse matrix A in CSR format
	cusparseCreateCsr(&matA, size, size, CSR_value_size,
		count_d, col_d,value_d,
		CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
		CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);
		// Create dense vector X
	cusparseCreateDnVec(&vecX,size, ans_CSR_d, CUDA_R_64F);
		// Create dense vector y
	cusparseCreateDnVec(&vecY, size, F_d, CUDA_R_64F);
		// allocate an external buffer if needed
	cusparseSpMV_bufferSize(
		handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
		&alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
		CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize);
	cudaMalloc(&dBuffer, bufferSize);

		// execute SpMV
	cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
		&alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
		CUSPARSE_SPMV_ALG_DEFAULT, dBuffer);

	cudaMemcpy(F,F_d, sizeof(double) * size, cudaMemcpyDeviceToHost);//从设备端取回结果
		// destroy matrix/vector descriptors
	cusparseDestroySpMat(matA);
	cusparseDestroyDnVec(vecX);
	cusparseDestroyDnVec(vecY);
	cusparseDestroy(handle);

	cudaFree(value_d);
	cudaFree(ans_CSR_d);
	cudaFree(count_d);
	cudaFree(col_d);
	cudaFree(F_d);
}