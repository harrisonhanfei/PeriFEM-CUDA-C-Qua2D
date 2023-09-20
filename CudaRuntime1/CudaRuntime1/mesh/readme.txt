左侧有缺口的版

每步步长 0.005
上边界竖向位移0.25
下边界竖向位移-0.25

计算信息
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