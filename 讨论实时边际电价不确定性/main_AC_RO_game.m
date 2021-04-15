%离散场景分布鲁棒优化，主函数，交流潮流，一次函数近似
%第九个时间断面开始，08:00-次日08:00
%列与约束生成算法
clear
clc
load DLMP_data
price_RT=1.5*[b_data(9:24);b_data(1:8)];
%电动汽车聚类数据(08:00-09:00为时段1)
EVdata=[6,40,15,10,24,105;6,32,16,2,9,132;3,24,8,11,23,89;6,24,12,13,22,97;6,40,25,1,8,101;6,40,16,12,23,175;3,24,8,10,24,35;10,40,20,2,8,88;10,40,18,11,24,112;10,64,25,11,23,66];
%不同类型电动汽车的比例
ratio=EVdata(:,6)/1000;
data_SP.ratio=ratio;%初始场景，即经验分布
for iter=1:2
    result_MP(iter)=RO_MP_AC_uncertain_game(data_SP,price_RT);%求解鲁棒主问题
    result_MP.obj%目标函数上界
    data_MP=result_MP(iter);%向子问题传递数据
    result_SP(iter)=RO_SP_AC_uncertain_game(data_MP,price_RT);%求解鲁棒子问题
    data_SP.ratio=[data_SP.ratio,result_SP(iter).ratio];%添加场景
    result_SP.obj%目标函数下界
    abs(result_MP(iter).obj-result_SP(iter).obj)/result_MP(iter).obj
    if abs(result_MP(iter).obj-result_SP(iter).obj)/result_MP(iter).obj<=1e-6%满足精度范围则退出
        break
    end
end
Pb_DA=data_MP.Pb_DA;
Pb_RT=data_MP.Pb_RT;
Ps_RT=data_MP.Ps_RT;
Pb_DA=[Pb_DA(17:24,1);Pb_DA(1:16,1)];
Pb_RT=[Pb_RT(17:24,1);Pb_RT(1:16,1)];
Ps_RT=[Ps_RT(17:24,1);Ps_RT(1:16,1)];
plot(Pb_DA)
hold on
plot(Pb_RT,'r')
plot(-Ps_RT,'g')
legend('日前购电','实时购电','实时售电')
xlim([1,24])
ylim([-550,550])