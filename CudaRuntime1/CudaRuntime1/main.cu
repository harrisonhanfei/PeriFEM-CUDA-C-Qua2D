//@author zjd
//2021.12.03   
//PeriFEM 二维程序示例v1.0  2021.12.27
#include"function.cuh"

int main()
{
	clock_t start, end, end1, end2, end3,end4, end5, end6;
	start = clock();//计时节点
	time_t curent_time,end_time;
	time(&curent_time);
	printf("GPU计算开始！！！\n");
	printf("@author zjd  ^·^\n");
	printf("当前时间：%s", ctime(&curent_time));
	printf("===================================\n");


	double*node = (double*)malloc(3 * NODE * sizeof(double));
	int*element = (int*)malloc(5 * ELE * sizeof(int));
	host_read_FEM_information(node, element);//读取FEM数据
	printf("===================================\n");
	int*relevant_element = (int*)malloc(ELE * skip_size * sizeof(int));//存储每个单元的相关单元
	
	float time_cal_rela_ele[1];
	printf("正在创建PeriFEM单元信息\r");
	host_find_relevant_element(node, element, relevant_element, time_cal_rela_ele);//寻找相关单元，即创建PE数据
	
	int big = 0, small = ELE, PeriFEM_num = 0;
	int big_half = 0, small_half = ELE, PeriFEM_num_half = 0;
	for (int i = 0; i < ELE; i++)
	{
		if (big < relevant_element[i*skip_size])big = relevant_element[i*skip_size];
		if (small > relevant_element[i*skip_size])small = relevant_element[i*skip_size];
		PeriFEM_num += relevant_element[i*skip_size];
		if (big_half < relevant_element[i*skip_size+1])big_half = relevant_element[i*skip_size+1];
		if (small_half > relevant_element[i*skip_size+1])small_half = relevant_element[i*skip_size+1];
		PeriFEM_num_half += relevant_element[i*skip_size+1];
	}
	printf("PeriFEM_ele_num :     %d\tPeriFEM_ele_num_half : %d\n", PeriFEM_num, PeriFEM_num_half);
	printf("big skip_size:        %d\tsmall skip_size:       %d\n", big, small);
	printf("big skip_size_half:   %d\tsmall skip_size_half:  %d\n", big_half, small_half);
	printf("===================================\n");

	int*node_to_node = (int*)malloc(NODE * node_to_node_skip * sizeof(int));//结点与节点的关系
	zeros(node_to_node_skip * NODE, node_to_node);

	int size[1];
	host_cal_node_to_node_relation(element, relevant_element, node_to_node, size);//计算节点间相关关系，为形成CSR指示矩阵条件
	int CSR_value_size = size[0];
	printf("CSR存储总纲需要value矩阵的大小为：%d\n", CSR_value_size);

	double*k_all_CSR_value = (double*)malloc(CSR_value_size * sizeof(double));//k_all的value
	int*k_all_CSR_count = (int*)malloc((NODE * 2 + 1) * sizeof(int));//行索引
	int*k_all_CSR_col = (int*)malloc(CSR_value_size * sizeof(int));//列索引

	int*PeriFEM_bond = (int*)malloc(ELE * skip_size * 17 * sizeof(int));//17个为一PeriFEM单元的键状态，第一个数记录键个数
	host_bond_initialize(ELE * skip_size * 17, PeriFEM_bond);//初始化键状态
	int*PeriFEM_bond_backups = (int*)malloc(ELE * skip_size * 17 * sizeof(int));//17个为一PeriFEM单元的键状态，第一个数记录键个数，备份键状态

	host_cal_CSR_count_col(CSR_value_size, k_all_CSR_col, k_all_CSR_count, node_to_node);//形成CSR指示矩阵
	free(node_to_node);

	double*k_PeriFEM_stiffness = (double*)malloc(ELE * skip_size_half * 136 * sizeof(double));
	
	float time_cal_k[1];
	host_calculate_PeriFEM_stiffness(node, element, relevant_element, k_PeriFEM_stiffness,time_cal_k);

	zeros_d(CSR_value_size, k_all_CSR_value);
	end1 = clock();
	host_k_all(element, relevant_element, k_all_CSR_count, k_all_CSR_col, k_all_CSR_value, k_PeriFEM_stiffness);//组总纲慢，因为与FEM比需要多组个相关单元倍数
	free(k_PeriFEM_stiffness);end2 = clock();
	printf("组装总纲用时：%.3fs。\n", (float)(end2 - end1) / CLOCKS_PER_SEC);
	
	printf("===================================\n");
	printf("开始迭代求解\n");
	printf("===================================\n");
	double*R = (double*)malloc(NODE * 2 * sizeof(double));
	double*ans_CSR = (double*)malloc(NODE * 2 * sizeof(double));
	
	float  time_bc = 0.0, time_CG = 0.0, time_bond = 0.0, time_update_kall = 0.0;

	int*PeriFEM_ele_broken = (int*)malloc(ELE * skip_size * sizeof(int));//存储每步的断键单元，每次最多坏ELE个，数组第一个表示主单元号，第二个计数，后面的表示从单元号
	double*k_all_CSR_value_ans = (double*)malloc(CSR_value_size * sizeof(double));//k_all的value备份
	double displacement = 0.0;
	FILE*print = fopen(".\\result\\计算中信息输出.txt", "a");
	while (true)
	{
		int inner = 0;
		displacement += 0.005;
		printf("当前位移步：%.4lf\n", displacement);
		fprintf(print, "当前位移步：%.4lf\n", displacement);
		while (true)
		{
			inner++;
			zeros_d(NODE * 2, R);
			zeros_d(NODE * 2, ans_CSR);
			
			copy(k_all_CSR_value_ans,k_all_CSR_value,CSR_value_size);
			float time_deal_bc[1];
			deal_bc(k_all_CSR_value_ans, k_all_CSR_count, k_all_CSR_col,CSR_value_size, R, node, displacement, time_deal_bc);//边界条件
			time_bc += time_deal_bc[0];
			
			int CG_num[1];
			end3 = clock();
			CG_iterative(k_all_CSR_value_ans, k_all_CSR_count, k_all_CSR_col, CSR_value_size, ans_CSR, R,CG_num);end4 = clock();//迭代（CG迭代求解
			time_CG += (float)(end4 - end3) / CLOCKS_PER_SEC;

			copy_bond(PeriFEM_bond_backups, PeriFEM_bond, ELE * skip_size * 17);//PeriFEM_bond-->PeriFEM_bond_backups
			int broken_ele_num[1] = { 0 };
			float time_up_bond[1];
			update_bond(node, element, relevant_element, ans_CSR, PeriFEM_bond, PeriFEM_ele_broken, broken_ele_num,time_up_bond);
			printf("当前迭代步%d    PeriFEM单元损伤数：%d\n", inner, broken_ele_num[0]);//更新键状态
			time_bond += time_up_bond[0];
			fprintf(print, "迭代步 %d  PeriFEM单元损伤数：%d\n", inner, broken_ele_num[0]);

			end5 = clock(); 
			update_k_all(node,element, relevant_element,k_all_CSR_count, k_all_CSR_col,
				k_all_CSR_value,PeriFEM_bond_backups, PeriFEM_bond, PeriFEM_ele_broken, broken_ele_num);//更新总纲
			end6 = clock();
			time_update_kall += (float)(end6 - end5) / CLOCKS_PER_SEC;
			if (broken_ele_num[0] == 0)break;
		}
		double*F = (double*)malloc(NODE * 2 * sizeof(double));//反力
		host_cal_F(k_all_CSR_count, k_all_CSR_col, k_all_CSR_value, CSR_value_size, ans_CSR, F);
		export_result_load(node, element, relevant_element, PeriFEM_bond,ans_CSR,F);
		free(F);
		if (displacement >= 0.25)break;
	}
	fclose(print);
	export_result_end(node, element, relevant_element, PeriFEM_bond,ans_CSR);

	end = clock();//计时节点
	time(&end_time);
	printf("计算完成！共用时：%.3fs。\n", (float)(end - start) / CLOCKS_PER_SEC);


	FILE*fp1;
	fp1 = fopen(".\\result\\CUDA计算信息输出.txt", "w");
	fprintf(fp1, "CUDA计算过程中信息输出\n");
	fprintf(fp1, "GPU计算开始！！！\n");
	fprintf(fp1, "@author zjd  ^·^\n");
	fprintf(fp1, "开始时间：%s", ctime(&curent_time));
	fprintf(fp1, "结束时间：%s", ctime(&end_time));
	fprintf(fp1, "===================================\n");
	fprintf(fp1, "GPU寻找相关单元用时: %f ms\n", time_cal_rela_ele[0]);
	fprintf(fp1, "PeriFEM_ele_num : %d\tPeriFEM_ele_num_half : %d\n", PeriFEM_num, PeriFEM_num_half);
	fprintf(fp1, "big skip_size:         %d\tsmall skip_size:          %d\n", big, small);
	fprintf(fp1, "big skip_size_half:   %d\tsmall skip_size_half:     %d\n", big_half, small_half);
	fprintf(fp1, "===================================\n");
	fprintf(fp1, "CSR存储总纲需要value矩阵的大小为：%d\n", CSR_value_size);
	fprintf(fp1, "GPU计算单刚用时:%.3f ms\n", time_cal_k[0]);
	fprintf(fp1, "CPU组装总纲用时:%.3f s\n", (float)(end2 - end1) / CLOCKS_PER_SEC);
	fprintf(fp1, "===================================\n");
	fprintf(fp1, "GPU更新边界用时:%.3f ms。\n", time_bc);
	fprintf(fp1, "迭代用时：%.3fs。\n",time_CG);
	fprintf(fp1, "GPU更新键状态用时:%.3f ms。\n", time_bond);
	fprintf(fp1, "===================================\n");
	fprintf(fp1, "CPU更新总纲用时:%.3f s。\n", time_update_kall);
	fprintf(fp1, "===================================\n");
	fprintf(fp1, "计算完成！共用时：%.3fs。", (float)(end - start) / CLOCKS_PER_SEC);
	fclose(fp1);


	free(node);
	free(element);
	free(relevant_element);
	free(k_all_CSR_col);
	free(k_all_CSR_count);
	free(k_all_CSR_value);
	free(PeriFEM_bond);
	free(PeriFEM_bond_backups);
	free(PeriFEM_ele_broken);
	free(k_all_CSR_value_ans);
	free(R);
	free(ans_CSR);

	return 0;
}