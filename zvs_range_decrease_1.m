clc
clear
%clf

k = 1.5;
p = 0.4;
I = 3; %charging current (A)
IN = 120/(8*5e-6*100000);
d1 = linspace(0,1,100);
i = I/IN;
d2_1 = (1-k*(1-d1)) / 2;    %condition 2 (iL(t1) < 0) 三电平高上升沿小于0，Q4关Q3开
d2_2 = (k-1-k*d1)/(2*k);    %condition 2 (iL(t2) > 0) 两电平上升沿大于0 Q6关Q5开
d2_3 = (k*(1-d1)-1 - i/2)/2;  %condition 1 (iL(t2) < 0) 三电平高上升沿小于0 Q4关Q3开
d2_4 = ((k-(1-i/2))/k);         %condition 1 (iL(t4) < 0 即 iL(t1) > 0) 两电平上升沿大于0 (下降沿小于0) Q6关Q5开

%转换为d_delta坐标

d_d2_1 = d2_1 + d1/2;
d_d2_2 = d2_2 + d1/2;
d_d2_3 = d2_3 + d1/2;
d_d2_4 = d2_4 + d1/2;
d_d_boud = 1 - d1/2;    %cond2 & cond3 的边界

d_d_boud2 = d1/2;      %condition1 & condition2 的边界

fy_d_d2_1 = @(d1)(1 + (i/2) - k*(1-d1)) / 2 + d1/2;   %cond2 边界1的解析表示
fy_d_d2_2 = @(d1)(2*k - 2*k*d1 - 2 + i)/(4*k) + d1/2;   %cond2 边界2的解析表示
fy_d_d_boud = @(d1)1 - d1/2;                %边界3的解析式
fun_cond2_1 = @(d1)fy_d_d2_1(d1) - fy_d_d2_2(d1);   %解1，2的交点(iL(t1)与iL(t2))
fun_cond2_4 = @(d1)fy_d_d2_1(d1) - fy_d_d_boud(d1);   %解1，3的焦点 (iL(t1)与iL(t3))

%在zvs缩小的时候，1与3的解是有的 1与2的解 不准确 要分别用cond1与cond2的边界与1，2求解(这句话可能是我之前理解的错误 暂时不对20191219)
fy_d_d_boud2 = @(d1)d1/2;
fun_cond2_11 = @(d1)fy_d_d2_1(d1) - fy_d_d_boud2(d1);   %解cond1与cond2边界线与1的解 相当于把fun_cond2_1扩展成两个点
fun_cond2_12 = @(d1)fy_d_d2_2(d1) - fy_d_d_boud2(d1);

%cond2 boundary line. 5lines with Im, 4lines without Im
fill_x_1 = fsolve(fun_cond2_1,[0,1]);             %取1，3交点的x坐标（因为解出来两个值相等）
fill_y_1 = (1 + (i/2) -k*(1-fill_x_1(1))) / 2 + fill_x_1(1)/2;     %解出1，3交点的y坐标
%fill_x_11 = fsolve(fun_cond2_11,[0,1]);             %取1与cond1&cond2边界的交点x坐标
%fill_y_11 = fill_x_11(1)/2;                         %解出y坐标，用边界线或者1的方程都可以
%fill_x_12 = fsolve(fun_cond2_12,[0,1]);
%fill_y_12 = fill_x_12(1)/2;
fill_x_4 = fsolve(fun_cond2_4,[0,1]);
fill_y_4 = 1 - fill_x_4(1)/2;
fill_x_2 = 0;
fill_y_2 = (k-1-k*fill_x_2)/(2*k) + fill_x_2/2;
fill_x_3 = 0;
fill_y_3 = 1 - fill_x_3/2;

fill_x_5 = d2_4;
fill_y_5 = 0;
fill_x_6 = d2_4;
fill_y_6 = (k*(1-fill_x_6)-1 - i/2)/2 + fill_x_6/2;
fill_y_7 = 0;                                       % 考虑Im时，与x轴交点会变
fill_x_7 = (2*fill_y_7 + i/2 + 1 - k)/(1-k);

%% 观察用d1和p代替d2后的电流最大值与d1取值的趋势
%I_max = 2*d1 - 2*k*(d1-1) - 2*sqrt(-d1.*d1 - p + 1);
d2_small = 1/2 - sqrt(- d1.*d1 - p + 1)/2 - d1/2;
d2_big = 1/2 + sqrt(- d1.*d1 - p + 1)/2 - d1/2;
d2_cond1 = - (p + 2*d1.*(d1 - 1))./(4*d1-4);
for cnt_x = 1:1:100
    if ~(isreal(d2_small(cnt_x)))
        d2_small(cnt_x) = nan;
    elseif d2_small(cnt_x) < 0
        d2_small(cnt_x) = nan;
    end
    if ~(isreal(d2_big(cnt_x)))
        d2_big(cnt_x) = nan;
    end
    if d2_cond1(cnt_x)>0
        d2_cond1(cnt_x) = nan;
    end
end

%将d2_small 与 d2_cond1进行合并，形成光滑曲线
d2_show = d2_small;
for cnt_x = 1:1:length(d1)
    if (isnan(d2_show(cnt_x)))
        d2_show(cnt_x) = d2_cond1(cnt_x);
    end
end

%寻找合并后small与big第一个nan值，赋予重合点值
for cnt_x = 1:1:length(d1)
    if (isnan(d2_show(cnt_x)))
        d2_show(cnt_x-1) = (1 - sqrt(1-p)) / 2; %因为取第一个nan值赋值不好看 会出现一个凸起，所以取上一个值。
        d2_big(cnt_x-1) = (1 - sqrt(1-p)) / 2;
        comb_flag = cnt_x;
        break;
    end
end

% 将df dd的所有坐标都转化成dd
d_d2_small = d2_small + d1/2;
d_d2_big = d2_big + d1/2;
d_d2_cond1 = d2_cond1 + d1/2;
d_d2_show = d2_show + d1/2;

I_max_small = 2 * ( k*(1-d1) + (2*d1 + 2*d2_small -1) );
I_max_big = 2 * ( k*(1-d1) + (2*d1 + 2*d2_big -1) );
I_max_t2 = 2*(k*(2*d2_small + d1 - 1) + 1);
I_max_t1 = 2*(k*(1-d1)+(2*d2_small-1));

d1_cond1_1 = -(2^(1/2)*(p*(k - 1))^(1/2) - 2*k + 2)/(2*(k - 1));
d2_cond1_1 = -(p + 2*d1_cond1_1*(d1_cond1_1 - 1))/(4*d1_cond1_1 - 4);
d_d2_cond1_1 = d2_cond1_1 + d1_cond1_1/2;
%if d2_cond1_1 > 0
%    d2_cond1_1 = 0;
%    d1_cond1_1 = 
%d1_cond1_2 = (2*k + 2^(1/2)*(p*(k - 1))^(1/2) - 2)/(2*(k - 1));
%I_max_cond1_1 = 4*d1_cond1_1 - (4*(p + 2*d1_cond1_1*(d1_cond1_1 - 1)))/(4*d1_cond1_1 - 4) - 2*k*(d1_cond1_1 - 1) - 2;
%I_max_cond1_2 = 4*d1_cond1_2 - (4*(p + 2*d1_cond1_2*(d1_cond1_2 - 1)))/(4*d1_cond1_2 - 4) - 2*k*(d1_cond1_2 - 1) - 2;
%I_max_cond1 = 4*d1 - (4*(p + 2*d1.*(d1 - 1)))./(4*d1 - 4) - 2*k*(d1 - 1) - 2;
I_max_cond1 = 2*(2*d1 + 2*d2_cond1 + k*(1-d1) - 1);     % iL(t3) is the peak current
I_max2_cond1 = -2*(k*(1-d1) - 1);                       % iL(t1) there is a minus that not iL(t4)
I_max3_cond1 = 2*(2*d2_cond1 + k*(d1-1) + 1);           % iL(t2) the other ZVS condition.
d = ( 1 - sqrt(1-p) ) / 2;
I_max_sps(1,1:100) = 2*(2*d - 1 + k);

%%%%%%%%找出cond1中的最大电流值%%%%%%%%%%
for cnt_2 = 1:1:length(d1)
    I_max_real_cond1(cnt_2) = max(max(I_max_cond1(cnt_2),I_max2_cond1(cnt_2)), I_max3_cond1(cnt_2));
end


%将I_max_cond1 与 I_max_small进行合并
I_max_show = I_max_small;
for cnt_x = 1:1:length(d1)
    if (isnan(I_max_show(cnt_x)))
        I_max_show(cnt_x) = I_max_cond1(cnt_x);
    end
end

%寻找合并后电流small与big第一个nan值，赋予重合点的电流值
for cnt_x = 1:1:length(d1)
    if (isnan(I_max_show(cnt_x)))
        I_max_show(cnt_x) = 2 * ( k*(1-d1(cnt_x)) + (2*d1(cnt_x) + 2*d2_show(cnt_x) -1) );
        I_max_big(cnt_x) =  2 * ( k*(1-d1(cnt_x)) + (2*d1(cnt_x) + 2*d2_big(cnt_x) -1) );
        break;
    end
end
%{
D1 = ((-(p - 1)*(k^2 - 2*k + 2))^(1/2)*(k - 1))/(k^2 - 2*k + 2);
D2 = 1/2 - (- D1^2 - p + 1)^(1/2)/2 - D1/2;
    if D2 < 0
        D1 = -(sqrt(2)*sqrt(p*(k-1)) - 2*k + 2) / (2*(k - 1));
        D2 = - (p + 2*D1*(D1-1)) / (4*D1-4);
        ilt2 = 2*( 2*D1 + 2*D2 + k*(1-D1) - 1 );
        ilt4 = -2 * (k*(1-D1) - 1);
        if ilt2<ilt4
            D1 = -(sqrt(2*p - 4*k*p + 1) - 4*k + 3) / (2*(2*k - 1));
            D2 = -(p + 2*D1*(D1-1)) / (4*D1 - 4);
        else
            D11 = (D1 + (k-1)/k) / 2;
            D2 = - (p + 2*D11*(D11-1)) / (4*D11-4);
            D1 = D11;
        end
    end
d_D2 = D2 + D1/2;
%}


%% 画图
figure (44)
hold on
box on
set(gca,'FontSize', 24, 'Fontname', 'Times New Roman');
xlabel('\itD_{\rmf}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
%yyaxis left
ylabel('\itD_{\rm\delta}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
%%%%%%%%%%%%%%%%以下为mode1，2，3边界%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fill([0,1,1], [0,0,0.5],'c');   % fill cond1
%fill([0,0,1], [0,1,0.5],'m');   % fill cond2
%fill([0,1,1], [1,0.5,1],'g');   % fill cond3
%plot(d1,d_d_boud,'--m','linewidth',2);  % cond2&cond3 边界
%plot(d1,d_d_boud2,'--k','linewidth',2); % cond1&cond2 边界
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%以下为zvs范围%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h1 = fill([fill_x_1(1),fill_x_2,fill_x_3,fill_x_4(1)], [fill_y_1,fill_y_2,fill_y_3,fill_y_4],'-y');
%h2 = fill([fill_x_5,fill_x_6,fill_x_7], [fill_y_5,fill_y_6,fill_y_7],'-y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
d2_5 = d1*0;
d_d2_5 = d2_5 + d1/2;

%plot(d1,d_d2_show,'-','linewidth',3);
%plot(d1,d_d2_big,'-','linewidth',3); %合并后的1，2移相角组合取值和cond3的取值
%plot(D1,d_D2,'r*'); % 生成全局最优点

%yyaxis right
%ylabel('\iti_{\rmpeak}{\rm(A)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
%plot(d1,I_max_show,'--',d1,I_max_big,'--'); %合并后的电流和cond3电流
%plot(d1,I_max_show,'--*'); 
%plot(d1,I_max_real_cond1,'-');
%plot(d1,I_max2_cond1);
%plot(d1,I_max_cond1,'--b','linewidth',2);
%plot(d1,I_max2_cond1,'-.g','linewidth',2);
%plot(d1,I_max3_cond1,':c','linewidth',2);
%plot(d1,I_max_sps,'-r','linewidth',2);
%plot([0,1],[0,0]);

%plot(d1,d_d2_small,d1,d_d2_big);

%plot(d1,d_d2_cond1);

%plot(D1,d_D2,'r*');

%plot(d1_cond1_1,d_d2_cond1_1,'b*');
hold on
%yyaxis right
%plot(d1,I_max_cond1,d1,I_max_small,d1,I_max_big);

%% whole range optimal curve

cnt = 1;
for P = 0.01:0.001:1
    D1 = (k-1)*sqrt( (1-P)*(k*k-2*k+2) ) / (k*k-2*k+2);
    D2 = 0.5 - sqrt(-D1*D1 - P + 1) / 2 - D1/2;
    if D2 < 0
        D1 = -(sqrt(2)*sqrt(P*(k-1)) - 2*k + 2) / (2*(k - 1));
        D2 = - (P + 2*D1*(D1-1)) / (4*D1-4);
        ilt2 = 2*( 2*D1 + 2*D2 + k*(1-D1) - 1 );
        ilt4 = -2 * (k*(1-D1) - 1);
        
        %*******************************************                              
        %如果要得到完全没有控制的an的曲线 要从此往下开始注释
        %*******************************************
        
        if ilt2<ilt4
            D1 = -(sqrt(2*P - 4*k*P + 1) - 4*k + 3) / (2*(2*k - 1));
            D2 = -(P + 2*D1*(D1-1)) / (4*D1 - 4);
        %else
        %    D11 = (D1 + (k-1)/k) / 2;
        %    D2 = - (P + 2*D11*(D11-1)) / (4*D11-4);
        %    D1 = D11;
        %end
        end
        
        %******************************************************
        %以下段落为考虑最小充电电流情况时对最优曲线进行zvs上的优化部分程序
        %******************************************************
        
        %if P >= 2*(k-1)/ ((3*k-2)*(3*k-2))
        if D2 > (k*(1.0-D1)-1 - i/2.0)/2.0 %不能得到反向正电流
            D1 = -(i - 4*k + (i^2 - 8*P + 8*k*P)^(1/2) + 4)/(4*(k - 1));
            D2 = -(P + 2*D1*(D1-1)) / (4*D1 - 4);   % 功率曲线和iL(t1)的交点
            if D1 < ((k-(1-i/2))/k) %left boundary
                
                D1 = -(sqrt(2)*sqrt(P*(k-1)) - 2*k + 2) / (2*(k - 1));
                D2 = - (P + 2*D1*(D1-1)) / (4*D1-4);
                ilt2 = 2*( 2*D1 + 2*D2 + k*(1-D1) - 1 );
                ilt4 = -2 * (k*(1-D1) - 1);
           
                if ilt2<ilt4
                    D1 = -(sqrt(2*P - 4*k*P + 1) - 4*k + 3) / (2*(2*k - 1));
                    D2 = -(P + 2*D1*(D1-1)) / (4*D1 - 4);
                end
            end
        end
        
        %******************************************************
        %以上段落为考虑最小充电电流情况时对最优曲线进行zvs上的优化部分程序
        %******************************************************
        %}
    end
    d1_wave(cnt) = D1;
    d2_wave(cnt) = (D2 + d1_wave(cnt)/2) ;
    cnt = cnt + 1;
end
%d2_wave(1) = 0;
%yyaxis left
plot(d1_wave,d2_wave,'-r','linewidth',3);
%plot(d1,d2_show+d1/2,'-',d1,d2_big+d1/2,'-');
x_left = 10;
x_right = 30;
%plot(d1(x_left:x_right),d2_show(x_left:x_right)+d1(x_left:x_right)/2,'.-');   %为了展示小范围的曲线
%plot(d1*180,(d2_show+d1/2)*180,'-','Linewidth', 3);
%plot(d1*180,(d2_show+d1/2)*180,'*','Linewidth', 3);
%plot([0,0],[0,90],'--r','Linewidth',3);

%% figure which do not use Df and D_delta, but just d1 and d2
%只讨论在mode1中的电流情况

figure(55)
hold on
box on
xlim([0,1]);
set(gca,'FontSize', 24, 'Fontname', 'Times New Roman');
xlabel('\itD_{\rmf}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
yyaxis left
ylabel('\itD_{\rm\delta}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
%plot(d1,d2_show+d1/2,'-',d1,d2_big+d1/2,'-');  %显示全范围的df，ddelta取值
plot(d1,d_d2_cond1,'-', 'linewidth', 3);  %只显示cond1中的df，ddelta取值
dd1 = -(sqrt(2)*sqrt(p*(k-1)) - 2*k + 2) / (2*(k - 1)); % dd1为不考虑sub-peak电流的iL(t3)最小值
dd2 = - (p + 2*dd1*(dd1-1)) / (4*dd1-4);
plot(dd1,dd2+dd1/2,'*','MarkerSize', 15, 'LineWidth', 10);
plot([dd1,dd1],[0,1],'-.k','linewidth',1);
%plot(D1,D2,'r*');
yyaxis right
ylim([-1,2]);
ylabel('\iti_{\rmL}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
plot(d1,I_max_cond1,'--m','linewidth',3);
for cnt_x = 1:1:100
    if (isnan(I_max_cond1(cnt_x)))
        I_max2_cond1(cnt_x) = NaN;
        I_max_sps(cnt_x) = NaN;
    end
end
%plot(d1,I_max_show,'-c','linewidth',3);  %此处将I_max2_cond1扩展成为了cond1，2中的small部分
plot(d1,I_max2_cond1,'-.g','linewidth',3);
%plot(d1,I_max3_cond1,':m','linewidth',3);  %为最小的电流 不予考虑，只用做检测用。
plot(d1,I_max_sps,'-r','linewidth',3);

% I_max_cond1 = iL(t3) peak current, I_max2_cond1 = iL(t1) left condition
% of ZVS, I_max3_cond1 = iL(t2) right condition of ZVS
%plot([0,1],[0,0]);
%legend('d\_d2\_small','d\_d2\_big','point','scaleplate','iL(t3)','iL(t1)','iL(t2)','iL\_SPS','y=0')

