%_________________________________________________________________________________
%  Blood-Sucking Leech Optimizer                                                               
%                                                                                                     
%  Developed in MATLAB R2023b                                                                  
%                                                                                                     
%  programming: Jianfu Bai                                                          
%                                                                                                     
%  e-Mail: Jianfu.Bai@UGent.be, magd.abdelwahab@ugent.be                                                               
%  Soete Laboratory, Department of Electrical Energy, Metals, Mechanical Constructions, and Systems, 
%  Faculty of Engineering and Architecture, Ghent University, Belgium                                                           
%                                                                                                                                                                                                                                                              
%  paper: Jianfu Bai, Hung Nguyen-Xuan, Elena Atroshchenko, Gregor Kosec, Magd Abdel Wahab, Blood-Sucking Leech Optimizer (2024).  
%____________________________________________________________________________________
clear all
close all
clc

tic
%benchmarksType = 1 for 23 Classical benchmark functions
%benchmarksType = 2 for CEC 2017
%benchmarksType = 3 for CEC 2019
%benchmarksType = 4 for CEC 2020
%benchmarksType = 5 for CEC 2022
benchmarksType = 1;

if benchmarksType == 1
    maxFunc = 23;
elseif benchmarksType == 2
    maxFunc = 30;
elseif benchmarksType == 3
    maxFunc = 10;
elseif benchmarksType == 4
    maxFunc = 10;
elseif benchmarksType == 5
    maxFunc = 12;    
else
    exit;
end

SearchAgents_no = 30;
Max_iteration= 1000;

tt_Max=1;
runs = 30;
RUN_Best_score_T_one=zeros(runs,maxFunc);
for fn1=1:1
tt=1;
disp(fn1);
RUN_Best_score_T=zeros(runs,tt_Max);
RUN_Best_score_Best = zeros(tt_Max);
RUN_Best_Score_Mean = zeros(tt_Max);
RUN_Best_Score_std = zeros(tt_Max);
while tt<=tt_Max
for fn = fn1:fn1
    
    Function_name=strcat('F',num2str(fn));
    if benchmarksType == 1
        [lb,ub,dim,fobj]=Get_Functions_details(Function_name);
    elseif benchmarksType == 2
        if fn == 2
            continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
        end
        [lb,ub,dim,fobj]=CEC2017(Function_name);
    elseif benchmarksType == 3
        [lb,ub,dim,fobj]=CEC2019(Function_name);
    elseif benchmarksType == 4   
        [lb,ub,dim,fobj]=CEC2020(Function_name);
    elseif benchmarksType == 5   
        [lb,ub,dim,fobj]=CEC2022(Function_name);        
    end
    
    Best_score_T = zeros(runs,1);
    AvgConvCurve = zeros(Max_iteration,1);
    Convergence_curve=zeros(runs,Max_iteration);

    Leeches_x=zeros(1,Max_iteration*SearchAgents_no);
    Leeches_y=zeros(1,Max_iteration*SearchAgents_no);
    %Execute the BSLO algorithm
    for run=1:runs
        % [Best_score,Best_pos,cg_curve,sh_x,sh_y]=HHO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
        [Best_score,Best_pos,cg_curve]=IBSLO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
        Best_score_T(run) = Best_score; 
        RUN_Best_score_T(run,tt)=Best_score;
    end
    Best_score_Best = min(Best_score_T);
    Best_Score_Mean = mean(Best_score_T);
    Best_Score_std = std(Best_score_T); 
    %output
    % format long
    % display([Function_name, ' Best:  ', num2str(Best_score_Best), '     ', 'Mean:  ', num2str(Best_Score_Mean), '     ', 'Std. Deviation:  ', num2str(Best_Score_std)]);
    RUN_Best_score_Best(tt) =Best_score_Best;
    RUN_Best_Score_Mean(tt) = Best_Score_Mean;
    RUN_Best_Score_std(tt) = Best_Score_std;   
    % % 
%     %draw search history
%     func_plot_CEC2017(Function_name);
%     title('Search history','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     xlabel('x_1','FontSize',21,'FontName','Times New Roman');
%     ylabel('x_2','FontSize',21,'FontName','Times New Roman');
%     zlabel([Function_name,'( x_1 , x_2 )'])
%     grid off
% 
%     sz = 25;
%     scatter(sh_x,sh_y,sz,'k','filled');
%     hold on;
%     % sz = 25;
%     % scatter(sh_x2,sh_y2,sz,'b','filled');
%     % hold on;    
%     sz = 50;
%     scatter(Best_pos(1,1),Best_pos(1,2),sz,'r','filled');
% 
%     set(gca,'xtick',lb:(ub-lb)/2:ub,'FontName','Times New Roman','FontSize',21);
%     set(gca,'ytick',lb:(ub-lb)/2:ub,'FontName','Times New Roman','FontSize',21);
%     % set(gca,'xtick',-5:5:10,'FontName','Times New Roman','FontSize',21);
%     % set(gca,'ytick',0:5:15,'FontName','Times New Roman','FontSize',21);
% figure;
end
tt=tt+1;
end
[Min_RUN_Best_Score_Mean,index]=min(RUN_Best_Score_Mean);
disp(index);
    Min_RUN_Best_score_Best =RUN_Best_score_Best(index(1));
    Min_RUN_Best_Score_Mean = RUN_Best_Score_Mean(index(1));
    Min_RUN_Best_Score_std = RUN_Best_Score_std(index(1)); 
 RUN_Best_score_T_one(:,fn1)=RUN_Best_score_T(:,index(1));
    format long
     display([Function_name, ' Best:  ', num2str(Min_RUN_Best_score_Best), '     ', 'Mean:  ', num2str(Min_RUN_Best_Score_Mean), '     ', 'Std. Deviation:  ', num2str(Min_RUN_Best_Score_std)]);

end
toc;
% 写入Excel文件
filename = 'IBSLO_23_results_500iter.xlsx';  % 定义Excel文件名
% 检查文件是否存在
if exist(filename, 'file')
    % 如果文件存在，则删除文件
    delete(filename);
    disp(['Deleted existing file: ', filename]);
else
    % 如果文件不存在，显示消息
    disp(['File not found, nothing to delete: ', filename]);
end
writematrix(RUN_Best_score_T_one, filename);