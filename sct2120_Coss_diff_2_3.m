clc
clear

a = 0.9636;
b = -0.3577;
c = -0.04179;
V3L = 360;
V2L = 240;
L = 125e3*9;    %这样设定使得电压电流单位为V和A，时间单位为ns
V_up = 360;

nums = xlsread('sct2120_Coss.xlsx');
vds = nums(:,1);
coss = nums(:,2).*(1e-3);    %由pf转化为nf
vds_inter = linspace(0.11,V_up,1000);    %插值的平均间隔横轴电压，要从1开始，不然前面几位会出现nan
coss_inter = interp1(vds,coss,vds_inter,'linear');
coss_2 = flip(coss_inter);
coss_il = coss_inter + coss_2;

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

%Edc_V2 = V1*qoss_il(end) - V2*qoss(end);
Edc = (2*V2L - V3L) * qoss(end);

il_0 = 0.1042;
il_end = sqrt(il_0^2 - 2*Edc/L);
t_freewheeling = il_end*L/(V2L);
vds8_0 = 0.00000001;
c8_vds8 = a*vds8_0^b + c;
c7_vds8 = a*(V3L - vds8_0)^b + c;
dvds8_0 = il_0/(c8_vds8-c7_vds8);

y0 = [vds8_0;dvds8_0];
[t,x] = ode45('vdp2_3',[0,2500],y0);
y = x(:,1);
dy = x(:,2);
figure(502)
hold on
box on
set(gca,'FontSize', 24, 'Fontname', 'Times New Roman');
xlabel('\itt {\rm(ns)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
ylabel('\itv_{\rmds} {\rm(V)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
plot(t,y,'linewidth',0.5);
for cnt = 1:length(y)
    if ~isreal(y(cnt))
    %if y(cnt) >= V1-0.999999999
        y(cnt) = NaN;
    end
end
%cnt = cnt-1;
%}
[m,p] = max(y);
t_left = t(p);
t_right = t_left + t_freewheeling;


il = (a*(y).^b+c + a*(V3L-y).^b+c).*dy;
%il = (a*(V3L - y).^b - a*y.^b).*dy;
figure(503)
hold on
box on
%plot(t,dy)
plot(t,il);