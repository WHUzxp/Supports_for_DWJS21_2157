function result_MP=RO_MP_DC_uncertain_game(data_SP,price_RT)
%% 基础数据
N=50;%电动汽车数量
%电价数据(08:00-09:00为时段1)
load DLMP_data
a_data=[a_data(9:24);a_data(1:8)];b_data=[b_data(9:24);b_data(1:8)];
%电动汽车聚类数据(08:00-09:00为时段1)
EVdata=[6,40,15,10,24,105;6,32,16,2,9,132;3,24,8,11,23,89;6,24,12,13,22,97;6,40,25,1,8,101;6,40,16,12,23,175;3,24,8,10,24,35;10,40,20,2,8,88;10,40,18,11,24,112;10,64,25,11,23,66];
utility=ones(10,24);
%不同类型电动汽车的比例
ratio=data_SP.ratio;%子问题得到的不同类型电动汽车的比例
[EVtype,Iterations]=size(ratio);
ratio_zhuanhuan=[];
for i=1:Iterations
    ratio_zhuanhuan=[ratio_zhuanhuan,ratio(:,i)];
end
%电动汽车状态矩阵
X=zeros(24,10);Y=zeros(24,10);
for i=1:10
    X(EVdata(i,4):EVdata(i,5),i)=1;%电动汽车停泊矩阵
end
%% 建模
price_EV=sdpvar(24,1);KKTsystem=[];
pch=sdpvar(24,10);%电动汽车充电
pdis=sdpvar(24,10);%电动汽车放电
S_EV=sdpvar(24,10);%电动汽车电量状态
for n=1:10
    obj_DR=price_EV'*(pch(:,n)-pdis(:,n))-utility(n,:)*(0.95*pch(:,n)-pdis(:,n)/0.95);%电动汽车目标函数
    C_EV=[0<=pch(:,n)<=ones(24,1)*EVdata(n,1).*X(:,n),0<=pdis(:,n)<=ones(24,1)*EVdata(n,1).*X(:,n),
        0.2*ones(24,1)*EVdata(n,2).*X(:,n)<=S_EV(:,n)<=0.95*ones(24,1)*EVdata(n,2).*X(:,n)];%电动汽车边界约束条件
    C_EV=[C_EV,S_EV(EVdata(n,4)+1:EVdata(n,5),n)==S_EV(EVdata(n,4):EVdata(n,5)-1,n)+0.95*pch(EVdata(n,4)+1:EVdata(n,5),n)-pdis(EVdata(n,4)+1:EVdata(n,5),n)/0.95,
        S_EV(EVdata(n,4),n)==EVdata(n,3)+0.95*pch(EVdata(n,4),n)-pdis(EVdata(n,4),n)/0.95,
        S_EV(EVdata(n,5),n)==0.95*EVdata(n,2)];%电动汽车电量约束条件
    ops=sdpsettings('kkt.dualbound',0);%不进行对偶边界估计
    [KKTsystem_single,details]=kkt(C_EV,obj_DR,price_EV,ops);%将电动汽车问题转化为KKT系统
    KKTsystem=[KKTsystem,KKTsystem_single];%合成KKT系统
    obj_EV(n)=utility(n,:)*(0.95*pch(:,n)-pdis(:,n)/0.95)-(details.b'*details.dual+details.f'*details.dualeq);
end
Pch=sdpvar(24,Iterations);%储能系统充电
Pdis=sdpvar(24,Iterations);%储能系统放电
S_ESS=sdpvar(24,Iterations);%储能系统电量状态
Pb_DA=sdpvar(24,1);%日前购电合同
price_DA=a_data.*Pb_DA+b_data;%日前DLMP
Pb_RT=sdpvar(24,Iterations);%实时购电量
Ps_RT=sdpvar(24,Iterations);%实时售电量
C_price=[mean(price_EV)<=mean(price_DA),0.8*price_DA<=price_EV<=1.2*price_DA];%零售电价约束条件
C_ESS=[0<=Pch<=250,0<=Pdis<=250,200<=S_ESS<=950,
    S_ESS(1)==500+0.95*Pch(1)-Pdis(1)/0.95,
    S_ESS(2:24)==S_ESS(1:23)+0.95*Pch(2:24)-Pdis(2:24)/0.95,
    S_ESS(24)==500];%储能系统约束条件
C_CS=[0<=Pb_DA<=500,0<=Pb_RT<=500,0<=Ps_RT<=500,Pb_DA*ones(1,Iterations)+Pb_RT+Pdis+N*pdis*ratio_zhuanhuan==Ps_RT+Pch+N*pch*ratio_zhuanhuan];%零售商约束条件
obj_CCG=sdpvar;%引入一个辅助变量用于CCG迭代
C_CCG=[];
for i=1:Iterations
    C_CCG=[C_CCG,obj_CCG<=sum(-(price_RT+0.001).*Pb_RT(:,i)+(price_RT-0.001).*Ps_RT(:,i))+N*obj_EV*ratio(:,i)];
end
%% 求解
ops=sdpsettings('solver','gurobi','gurobi.FeasibilityTol',1e-9,'gurobi.IntFeasTol',1e-9,'gurobi.MIPGap',1e-9,'gurobi.OptimalityTol',1e-9);%求解器参数,MILP问题
Constraints=[C_price,C_ESS,C_CS,KKTsystem,C_CCG];%总的约束条件
obj=-price_DA'*Pb_DA+obj_CCG;%总的目标函数（零售商的收益）
result=optimize(Constraints,-obj,ops)%求解最大化问题
price_EV=double(price_EV);pch=double(pch);pdis=double(pdis);S_EV=double(S_EV);Pch=double(Pch);Pdis=double(Pdis);S_ESS=double(S_ESS);Pb_DA=double(Pb_DA);Pb_RT=double(Pb_RT);Ps_RT=double(Ps_RT);
result_MP.obj=double(obj);result_MP.pch=pch;result_MP.pdis=pdis;result_MP.price_EV=price_EV;result_MP.Pb_DA=Pb_DA;price_DA=double(price_DA);
result_MP.Pb_RT=Pb_RT;result_MP.S_EV=S_EV;result_MP.Pch=Pch;result_MP.Pdis=Pdis;result_MP.S_ESS=S_ESS;result_MP.Ps_RT=Ps_RT;result_MP.price_DA=price_DA;
end
