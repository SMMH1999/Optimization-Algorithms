%  Sine Cosine Algorithm (SCA)  
%
%  Source codes demo version 1.0                                                                      
%                                                                                                     
%  Developed in MATLAB R2011b(7.13)                                                                   
%                                                                                                     
%  Author and programmer: Seyedali Mirjalili                                                          
%                                                                                                     
%         e-Mail: ali.mirjalili@gmail.com                                                             
%                 seyedali.mirjalili@griffithuni.edu.au                                               
%                                                                                                     
%       Homepage: http://www.alimirjalili.com                                                         
%                                                                                                     
%  Main paper:                                                                                        
%  S. Mirjalili, SCA: A Sine Cosine Algorithm for solving optimization problems
%  Knowledge-Based Systems, DOI: http://dx.doi.org/10.1016/j.knosys.2015.12.022

% This function draws the benchmark functions

function func_plot(func_name)

[lb,ub,dim,fobj]=Get_Functions_details(func_name);

switch func_name 
    case 'F1' 
        x=-100:2:100; y=x; %[-100,100]
        
    case 'F2' 
        x=-10:0.2:10; y=x; %[-10,10]
        
    case 'F3' 
        x=-100:2:100; y=x; %[-100,100]
        
    case 'F4' 
        x=-100:2:100; y=x; %[-100,100]
    case 'F5' 
        x=-30:0.6:30; y=x; %[-30,30]
    case 'F6' 
        x=-100:2:100; y=x; %[-100,100]
    case 'F7' 
        x=-1.28:0.0256:1.28;  y=x;  %[-1.28,1.28]
    case 'F8' 
        x=-500:10:500;y=x; %[-500,500]
    case 'F9' 
        x=-5.12:0.1024:5.12;   y=x; %[-5.12,5.12]    
    case 'F10' 
        x=-32:0.64:32; y=x;%[-32,32]
    case 'F11' 
        x=-600:12.5:600; y=x;%[-600,600]
    case 'F12' 
        x=-50:1:50; y=x;%[-50,50]
    case 'F13' 
        x=-50:1:50; y=x;%[-50,50]
    case 'F14' 
        x=-65.536:1.31072:65.536; y=x;%[-65.536,65.536]
    case 'F15' 
        x=-5:0.1:5; y=x;%[-5,5]
    case 'F16' 
        x=-5:0.1:5; y=x;%[-5,5]
    case 'F17' 
        x=-5:0.15:10; y=0:0.15:15;
    case 'F18' 
        x=-2:0.04:2; y=x;%[-2,2]
    case 'F19' 
        x=0:0.01:1; y=x;%[0,1]
    case 'F20' 
        x=0:0.01:1; y=x;%[0,1]        
    case 'F21' 
        x=0:0.1:10; y=x;%[0,10]
    case 'F22' 
        x=0:0.1:10; y=x;%[0,10]    
    case 'F23' 
        x=0:0.1:10; y=x;%[0,10] 
end    

    

L=length(x);
f=[];

for i=1:L
    for j=1:L
        if strcmp(func_name,'F15')==0 && strcmp(func_name,'F19')==0 && strcmp(func_name,'F20')==0 && strcmp(func_name,'F21')==0 && strcmp(func_name,'F22')==0 && strcmp(func_name,'F23')==0
            f(i,j)=fobj([x(i),y(j)]);
        end
        if strcmp(func_name,'F15')==1
            f(i,j)=fobj([x(i),y(j),0,0]);
        end
        if strcmp(func_name,'F19')==1
            f(i,j)=fobj([x(i),y(j),0]);
        end
        if strcmp(func_name,'F20')==1
            f(i,j)=fobj([x(i),y(j),0,0,0,0]);
        end       
        if strcmp(func_name,'F21')==1 || strcmp(func_name,'F22')==1 ||strcmp(func_name,'F23')==1
            f(i,j)=fobj([x(i),y(j),0,0]);
        end          
    end
end

surfc(x,y,f,'LineStyle','none');
title('Parameter space','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
xlabel('x_1','FontSize',21,'FontName','Times New Roman');
ylabel('x_2','FontSize',21,'FontName','Times New Roman');
zlabel(func_name,'FontSize',21,'FontName','Times New Roman');
set(gca,'xtick',lb:(ub-lb)/2:ub,'FontName','Times New Roman','FontSize',21);
set(gca,'ytick',lb:(ub-lb)/2:ub,'FontName','Times New Roman','FontSize',21);
% set(gca,'xtick',-5:5:10,'FontName','Times New Roman','FontSize',21);
% set(gca,'ytick',0:5:15,'FontName','Times New Roman','FontSize',21);
% set(gca,'ztick',0:60:120,'FontName','Times New Roman','FontSize',21);
% 绘制三维等高线图
figure;
contour(x,y,f);
% colorbar;
hold on;
end

