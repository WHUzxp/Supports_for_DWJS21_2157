%ÄâºÏÏµÊý
clear
clc
load AC_nolimit_DA
a_data=zeros(24,1);b_data=zeros(24,1);
x=10*[0:50];
for t=1:24
    a=sdpvar;b=sdpvar;
    y=DLMP_data(t,:);
    obj=sum((a*x+b-y).^2);
    optimize([],obj);
    a_data(t)=double(a);b_data(t)=double(b);
end