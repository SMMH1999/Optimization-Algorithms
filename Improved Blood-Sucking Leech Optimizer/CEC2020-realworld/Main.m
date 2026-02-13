%_________________________________________________________________________________
%  Blood-Sucking Leech optimizer: a new bio-inspired meta-heuristic optimization algorithm                                                               
%                                                                                                     
%  Developed in MATLAB R2022a                                                                  
%                                                                                                     
%  programming: Jianfu Bai                                                          
%                                                                                                     
%  e-Mail: Jianfu.Bai@UGent.be, magd.abdelwahab@ugent.be                                                               
%  Soete Laboratory, Department of Electrical Energy, Metals, Mechanical Constructions, and Systems, 
%  Faculty of Engineering and Architecture, Ghent University, Belgium                                                           
%                                                                                                                                                                                                                                                              
%  paper: Jianfu Bai, Samir Khatir, Laith Abualigah, Magd Abdel Wahab, Blood-Sucking Leech optimizer: a new bio-inspired meta-heuristic optimization algorithm (2023).  
%____________________________________________________________________________________

clear all
clc
format long
global initial_flag
%benchmarksType = 1 for spring
%benchmarksType = 2 for vessel
%benchmarksType = 3 for weldbeam
%benchmarksType = 4 for speed reducer
%benchmarksType = 5 for 3-bar truss
fn = 1;

tt_Max=1;
runs = 30;
RUN_Best_score_T_one=zeros(runs,33);
for fn=15:20
if fn==24
    tt_Max=2;
elseif fn==33
    tt_Max=2;
else
    tt_Max=10;
end
disp(['F',num2str(fn)]);
RUN_Best_score_T=zeros(runs,tt_Max);
RUN_Best_score_Best = zeros(tt_Max);
RUN_Best_Score_Mean = zeros(tt_Max);
RUN_Best_Score_std = zeros(tt_Max);

t=1;
while t<=tt_Max
initial_flag = 0;
Function_name=strcat('F',num2str(fn));
[lb,ub,dim,fobj]=Get_cec20_func(Function_name); 
Max_iteration= 1000;
   if dim <= 10
       MaxFES = 1e5;       
   elseif dim <= 30 ||dim > 10
       MaxFES = 2e5;
   elseif dim <= 50 ||dim > 30
       MaxFES = 4e5;
   elseif dim <= 150||dim > 50
       MaxFES = 8e5;
   else
       MaxFES = 1e6;
   end
SearchAgents_no=MaxFES/Max_iteration;   
% Calling algorithm
Best_score_T = zeros(runs,1);
AvgConvCurve = zeros(Max_iteration,1);
Convergence_curve=zeros(runs,Max_iteration);
Best_position=zeros(runs,dim);
%display (['Function:   ', num2str(fn)]);
for run=1:runs
    %rng('shuffle');
    [Best_score,Best_pos,cg_curve]=AGJO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
    Best_score_T(run) = Best_score;
%     Convergence_curve(run,:)=cg_curve;
    Best_position(run,:)=Best_pos;
    RUN_Best_score_T(run,t)=Best_score;
end
    Best_score_Best = min(Best_score_T);
    Best_Score_Mean = mean(Best_score_T);
    Best_Score_std = std(Best_score_T); 
    %output
    % format long
    % display([Function_name, ' Best:  ', num2str(Best_score_Best), '     ', 'Mean:  ', num2str(Best_Score_Mean), '     ', 'Std. Deviation:  ', num2str(Best_Score_std)]);
    RUN_Best_score_Best(t) =Best_score_Best;
    RUN_Best_Score_Mean(t) = Best_Score_Mean;
    RUN_Best_Score_std(t) = Best_Score_std;  

t=t+1;
end
[Min_RUN_Best_Score_Mean,index]=min(RUN_Best_Score_Mean);
disp(index);
    Min_RUN_Best_score_Best =RUN_Best_score_Best(index(1));
    Min_RUN_Best_Score_Mean = RUN_Best_Score_Mean(index(1));
    Min_RUN_Best_Score_std = RUN_Best_Score_std(index(1)); 
 RUN_Best_score_T_one(:,fn)=RUN_Best_score_T(:,index(1));
    format long
     display([Function_name, ' Best:  ', num2str(Min_RUN_Best_score_Best), '     ', 'Mean:  ', num2str(Min_RUN_Best_Score_Mean), '     ', 'Std. Deviation:  ', num2str(Min_RUN_Best_Score_std)]);

end
% 写入Excel文件
filename = 'AGJO_RC_CEC2020_results.xlsx';  % 定义Excel文件名
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