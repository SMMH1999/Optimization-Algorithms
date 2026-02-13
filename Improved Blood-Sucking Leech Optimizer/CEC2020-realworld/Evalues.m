clear all
close all
clc
Ev=zeros(1,500);
t=1;   
Max_iter=1000;
CM=chaos(9,Max_iter,1); 
while t<1001
aa=rand();
     if aa<0.5
     a=8-1*(-(t/Max_iter)^2+1);
     else
     a=8-7*(-(t/Max_iter)^2+1);
     end
    TD=a*(1-(t/Max_iter));
    r1=2*rand()-1;
    E=TD*r1*(CM(t)-(t/Max_iter));
    Ev(t)=E;
    t=t+1;
end
plot(Ev,'red');
xlabel('Iteration','FontSize',18,'FontName','Times New Roman');
ylabel('Perceived Distance (PD)','FontSize',18,'FontName','Times New Roman');
ax = gca; % 获取当前坐标轴对象
ax.XAxis.FontSize = 18; % 设置X轴标签字体大小为12
ax.XAxis.FontName = 'Times New Roman'; % 设置X轴标签字体为Arial
ax.YAxis.FontSize = 18; % 设置Y轴标签字体大小为12
ax.YAxis.FontName = 'Times New Roman'; % 
% hold on;
% t=1;   
% while t<1001
% aa=rand();
%      if aa<0.5
%      a=8-1*(-(t/Max_iter)^2+1);
%      else
%      a=8-7*(-(t/Max_iter)^2+1);
%      end
%     TD=a*(1-(t/Max_iter));
%     r1=2*rand()-1;
%     E=TD*r1*CM(t);
%     Ev(t)=E;
%     t=t+1;
% end
% plot(Ev,'g');