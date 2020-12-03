clc
clear

V1 = 240;
V2 = 360;
L = 125e-6;
V_up = 700;

nums = xlsread('sct2120_Eoss.xlsx');
vds = nums(:,1);
eoss = nums(:,2);
vds_inter = linspace(0,700,10000);
eoss_inter = interp1(vds,eoss,vds_inter,'linear');
figure(401)
hold on
box on
xlim([0,V_up]);
%ylim([-1100,1100]);
set(gca,'FontSize', 24, 'Fontname', 'Times New Roman');
xlabel('\itv_{\rmds}{\rm(V)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
ylabel('\itE_{\rmoss} {\rm(\itv)}{\rm(uJ)}', 'FontSize', 26.4, 'Fontname', 'Times New Roman');
%plot(vds,eoss,vds_inter,eoss_inter);
plot(vds, eoss,'linewidth',3);
plot(vds_inter, eoss_inter,'linewidth',3);

%% another method for verification
nums2 = xlsread('sct2120_Coss_1000V.xlsx');
vds2 = nums2(:,1);
coss2 = nums2(:,2).*(1e-3); %pf->nf
vds_inter2 = linspace(0.11,V_up,1000);
dv = vds_inter2(2) - vds_inter(1);
coss_inter2 = interp1(vds2, coss2, vds_inter2, 'linear');
eoss2(1) = vds_inter2(1)*coss_inter2(1)*dv;
for x = 2:1:length(coss_inter2)
    eoss2(x) = vds_inter2(x)*coss_inter2(x)*dv+eoss2(x-1);
end
plot(vds_inter2, eoss2*1e-3,'linewidth',3);



qoss(1) = coss_inter2(1)*dv;
for x = 2:1:length(coss_inter2)
    qoss(x) = coss_inter2(x)*dv + qoss(x-1);
end

eoss1(1) = vds_inter2(1) * qoss(1);
for x = 2:1:length(qoss)
    eoss1(x) = vds_inter2(x) * (qoss(x)-qoss(x-1)) + eoss1(x-1);
end
plot(vds_inter2, eoss1*1e-3,'linewidth',3);

eoss3(1) = vds_inter2(1) * (qoss(2)-qoss(1));
for x = 2:1:(length(qoss)-1)
    eoss3(x) = vds_inter2(x) * (qoss(x+1) - qoss(x)) + eoss3(x-1);
end
l = length(qoss);
eoss3(l) = vds_inter2(l)*(qoss(l) - qoss(l-1)) +eoss3(l-1);
plot(vds_inter2, eoss3*1e-3,'linewidth',3);
%}