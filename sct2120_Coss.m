clc
clear

V1 = 240;
V2 = 360;
L = 125e-6;
V_up = 700;

nums = xlsread('sct2120_Coss_1000V.xlsx');
vds = nums(:,1);
coss = nums(:,2).*(1e-3);    %由pf转化为nf
vds_inter = linspace(0.11,V_up,1000);    %插值的平均间隔横轴电压，要从1开始，不然前面几位会出现nan
coss_inter = interp1(vds,coss,vds_inter,'linear');
coss_2 = flip(coss_inter);
coss_il = coss_inter + coss_2;

%%%%%%%%%%%%%%%%对Coss拟合%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 3;
%p = polyfit(vds_inter,coss_inter,n);
p = polyfit(vds_inter,log(coss_inter),1);

coss_fit = exp(polyval(p,vds_inter));
%{
coss_fit = p(n+1);
for cnt = 1:1:n
    coss_fit = coss_fit + p(cnt) * vds_inter.^(n-cnt+1);
end
%}
%coss_fit = p(1)*vds_inter.^4 + p(2)*vds_inter.^3 + p(3)*vds_inter.^2 ...,
%            +p(4)*vds_inter.^1 + p(5);
%a = 1.06;
%b = -0.3603;
%c = -0.0798;
a = 0.9636;
b = -0.3577;
c = -0.04179;
coss_fit = a * vds_inter.^b + c;

%%%%%%%%%%%%%用来检测插值后的曲线和翻转曲线%%%%%%%%%%%%%%%%%

figure(301)
hold on
box on
xlim([0,V_up]);
set(gca,'FontSize', 24, 'Fontname', 'Times New Roman');
xlabel('\itv_{\rmds} {\rm(V)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
ylabel('\itC_{\rmoss} {\rm(nF)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
plot(vds,coss, 'linewidth',3);
plot(vds_inter, coss_inter,'--','linewidth',3);
plot(vds_inter, coss_fit,'-','linewidth',3);
%plot(vds_inter, coss_fit,'--','linewidth',3);
plot(vds_inter, coss_2, 'linewidth',3);
%plot(vds_inter, coss_il, 'linewidth',3);
%legend('{\itC}_3({\itv}) from datasheet', '{\itC}_3({\itv}) by interploation',...
%    '{\itC}_4({\itv})', '{\itC}_3({\itv}) + {\itC}_4({\itv})');
%}

l = length(vds_inter);

delta_q(1) = coss_inter(1) * (vds_inter(1)-0);
delta_q_il(1) = coss_il(1) * (vds_inter(1) - 0);

for cnt = 1:1:l-1
    delta_q(cnt) = coss_inter(cnt) * (vds_inter(cnt+1) - vds_inter(cnt));
    delta_q_il(cnt) = coss_il(cnt) * (vds_inter(cnt+1) - vds_inter(cnt));
end
delta_q(l) = coss_inter(l)*(vds_inter(2) - vds_inter(1));
delta_q_il(l) = coss_il(l)*(vds_inter(2) - vds_inter(1));

qoss(1) = delta_q(1);
qoss_il(1) = delta_q_il(1);
for cnt = 2:1:l
    qoss(cnt) = qoss(cnt-1) + delta_q(cnt);
    qoss_il(cnt) = qoss_il(cnt-1) + delta_q_il(cnt);
end
dv = vds_inter(2) - vds_inter(1);
qoss2(1) = coss_2(1) * dv;
for x = 2:1:length(coss_2)
    qoss2(x) = coss_2(x)*dv + qoss2(x-1);
end


x = 0;
for il_0 = 0:0.001:0.5
    x = x + 1;
    fx_sum(x) = 0;
    il_0_wave(x) = il_0;
    for cnt = 2:1:l
        il_v(cnt) = sqrt(il_0^2 + (2/L)*(qoss_il(cnt)*V2 - qoss(cnt)*V1)*1e-9); 
        fx = coss_il(cnt)/il_v(cnt)*(vds_inter(cnt) - vds_inter(cnt-1));
        fx_sum(x) = fx_sum(x) + fx;
    end
end
Edc = V2*qoss_il - V1*qoss;

Eoss(1) = vds_inter(1) * qoss(1);
for x = 2:1:length(qoss)
    Eoss(x) = vds_inter(x) * (qoss(x) - qoss(x-1)) + Eoss(x-1);
end

Eoss2(1) = vds_inter(1) * qoss2(1);
for x = 2:1:length(qoss2)
    Eoss2(x) = vds_inter(x) * (qoss2(x) - qoss2(x-1)) + Eoss2(x-1);
end

Eoss_il(1) = vds_inter(1) * qoss_il(1);
for x = 2:1:length(qoss_il)
    Eoss_il(x) = vds_inter(x) * (qoss_il(x) - qoss_il(x-1)) + Eoss_il(x-1);
end


%Eoss = vds_inter .* qoss;
figure(401)
hold on
box on
xlim([0,V_up]);
%ylim([-1100,1100]);
set(gca,'FontSize', 24, 'Fontname', 'Times New Roman');
xlabel('\itv_{\rmds}{\rm(V)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
ylabel('\itE_{\rmoss} {\rm(\itv)}{\rm(uJ)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
x = 0;
plot(vds_inter, Eoss*1e-3,'linewidth',3);
plot(vds_inter, Eoss2*1e-3,'linewidth',3);
plot(vds_inter, Eoss_il*1e-3,'linewidth',3);
plot(vds_inter, (Eoss+Eoss2)*1e-3,'--','linewidth',3);
legend('{\itE}_{oss1}','{\itE}_{qoss2}','{\itE}_{oss\_il}','{\itE}_{oss1} + {\itE}_{qoss2}')

figure(300)
hold on 
box on
xlim([0,V_up]);
set(gca,'FontSize', 24, 'Fontname', 'Times New Roman');
xlabel('\itv_{\rmds} {\rm(V)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
ylabel('\itQ_{\rmoss} {\rm(nC)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
plot(vds_inter,qoss, 'linewidth',3);
plot(vds_inter,qoss_il,'linewidth',3);
plot(vds_inter,qoss2,'linewidth',3);
%plot(vds_inter,il_v);
%plot(vds_inter,Edc);
legend('{\itQ}_{oss}({\itv})','{\itQ}_{oss\_iL}({\itv})','{\itQ}_{oss2}({\itv})')

figure(302)
hold on
box on
xlim([0,0.5]);
set(gca,'FontSize', 24, 'Fontname', 'Times New Roman');
xlabel('\iti_{\rmL}{\rm(0)} {\rm(A)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
ylabel('\itt_{\rmd} {\rm(ns)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
plot(il_0_wave, fx_sum, 'linewidth', 3);
%plot(vds,coss_il);

figure(303)
hold on
box on
xlim([0,V_up]);
%ylim([-1100,1100]);
set(gca,'FontSize', 24, 'Fontname', 'Times New Roman');
xlabel('\itv_{\rmds}{\rm(V)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
ylabel('\itE_{\rmdc} {\rm(\itv)}{\rm(nJ)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
x = 0;

%e_dc =  V_up*qoss - 120*qoss_il;   %k>1
e_dc =  120 * qoss_il - 120 * qoss;  %k<1

plot(vds_inter, e_dc,'linewidth',3);

plot([0,V_up],[0,0],'--k')



%}