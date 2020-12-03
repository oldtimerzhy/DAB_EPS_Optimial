# DAB_EPS_Optimial
MATLAB
zvs_range: 生成fig2,3,4， 给定k，p的软开关范围和取到最优质值。fig2为cond2，fig3为cond1，fig4为合并起来的软开关范围。
zvs_range_figure: 给定k或p，生成全范围的移相角取值和电流大小。 同时还有一种输出是以角度（乘180度）表示，将small与cond1合并。颜色RBG: 1: FFFFCC, 2: CCFFFF, 3: FFCCCC
traversal_zvs: 遍历程序，用来遍历是否cond1或cond2的所有移相角取值在给定范围内满足软开关。
optimal_Imax_cond1_cond2: 生成最优取值下的Imax电流峰值3D图形。拟加入rms值的3D图
zvs_range_cond1:专门负责研究cond1各个电流与sps峰值电流和eps自身rms电流的大小关系，选择最佳值。
Imax_SPS: 为了展示随着k的增大，SPS的峰值电流会增大
process_equation: 主要是为cond1的公式求解
single_test: 是在实验中发现cond1中有时候并不能实现两个峰值相等，通过调整D1的大小来手动的降低其中一侧的电流peak值，使取值向左走。这个程序是单纯针对实验中110/120的情况，其他情况需要考虑p，k取值而导致计算公式上的不同。
Irms_figure: 用来画出D1和P给定情况下的rms电流情况
map_generation: 用来生成不同功率等级下的电流，同时还可以得到具体的电流值和功率值，可用来画map与实验数据对比。
whole_range_optimal_curve: 在全范围内画出mode1,2的软开关范围，然后根据全功率范围最优值进行取点划线。
zvs_range_decrease: 标出来每个边界对应的哪个点电流，与zvs_range没有代码上改动。
zvs_range_decrease_1 :在figure4中生成减小的zvs范围。同时还能观察cond1中的各个peak电流大小关系。！！！！目前计算优化曲线的程序！！！
calculated_verification: 用来生成与两篇参考文献相比较的peak电路曲线。
calculated_verification2: 用来生成与B.Zhao文章中的peak电流比较。
min_charging_current: 利用公式计算最小充电电流和充电时间。
cal_result_with_parastic_wave: 考虑寄生电阻情况下对eps和sps各个顶点电流的修正。
calculated_verification: 和an的论文进行对照。
static_error_component: 用于验证10kWDAB的静态误差，其中包含逐点的D1、D2移相角的数值获得。
minus_Im_control_1kW: 为了获得1kWDAB的负Im控制中，图像和逐点的D1、D2移相角数值。
cal_result_with_parastic_numerical: 用来得到考虑和不考虑寄生电阻下的ZVS右边界曲线，用于Im控制画图。
boundary_and_power: 用来画mode1，2，3的区域和边界，同时画与k无关的power范围。
Journal_zvs_range_decrease_1: 用来画杂志论文中，不同k值混合的优化点曲线。
calculated_verification_only_sps: 用来生成只有SPS时，功率从大到小，不同k下的peak current。
calculated_verification_peak_current_with_k: 用来生成在Df和Ddelta坐标系下的，不用k在相同p参考下的各个peak current.

Numerical Analysis parts:
vdp.m: 用来列出1，3象限的微分方程。
vdp_2.m: 用来列出4，2象限的微分方程。
vdp_2_3.m: 利用新的推导来列的4，2象限微分方程，其结果公式和vdp_2可推导成相同，其中V1，V2换成了V2L和V3L表示。
sct2120_Coss_diff.m: 数值分析用来求解vdp.m的微分方程。
sct2120_Coss_diff_2.m: 数值分析用来求解vdp_2.m的微分方程。
sct2120_Coss_diff_2_3.m: 数值分析用来求解vdp_2_3的微分方程，其中V1，V2换成了V2L和V3L表示。
