clear all
close all
clc
format long
global initial_flag
%benchmarksType = 1 for spring
%benchmarksType = 2 for vessel
%benchmarksType = 3 for weldbeam
%benchmarksType = 4 for speed reducer
%benchmarksType = 5 for 3-bar truss
fn = 1;

tt_Max=10;
runs = 1;

for fn=15:15
if fn==24
    tt_Max=5;
elseif fn==33
    tt_Max=5;
else
    tt_Max=10;
end
        if fn == 23
            continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
        end
        if fn == 25
            continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
        end
        if fn == 26
            continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
        end
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
    [Best_score,Best_pos,cg_curve]=HOA_v2(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
    Best_score_T(run) = Best_score;
%     Convergence_curve(run,:)=cg_curve;
    Best_position(run,:)=Best_pos;
end
    Best_score_Best = min(Best_score_T);
    Best_Score_Mean = mean(Best_score_T);
    Best_Score_std = std(Best_score_T); 
    %output
    format long
    display([Function_name, ' Best:  ', num2str(Best_score_Best), '     ', 'Mean:  ', num2str(Best_Score_Mean), '     ', 'Std. Deviation:  ', num2str(Best_Score_std)]);

    x=zeros(1,Max_iteration);
    for tt=1:1000
        x(tt)=tt;
    end
    % plot(x,AvgConvCurve,'Color','b','LineWidth',3);
    % title('Average fitness','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
    % xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
    % set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);

    figure;
    semilogy(x,cg_curve,'Color','b','LineWidth',3);
    title('Convergence curve','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
    xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
    set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
    t=t+1;
end

end


% tt_Max=10;
% runs = 1;
% 
% for fn=24:33
% if fn==24
%     tt_Max=10;
% elseif fn==33
%     tt_Max=10;
% else
%     tt_Max=5;
% end
%         if fn == 23
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 25
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 26 || fn == 27 || fn == 28 || fn == 29 || fn == 30 || fn == 31 || fn == 32
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
% t=1;
% while t<=tt_Max
% initial_flag = 0;
% Function_name=strcat('F',num2str(fn));
% [lb,ub,dim,fobj]=Get_cec20_func(Function_name); 
% Max_iteration= 1000;
%    if dim <= 10
%        MaxFES = 1e5;       
%    elseif dim <= 30 ||dim > 10
%        MaxFES = 2e5;
%    elseif dim <= 50 ||dim > 30
%        MaxFES = 4e5;
%    elseif dim <= 150||dim > 50
%        MaxFES = 8e5;
%    else
%        MaxFES = 1e6;
%    end
% SearchAgents_no=MaxFES/Max_iteration;   
% % Calling algorithm
% Best_score_T = zeros(runs,1);
% AvgConvCurve = zeros(Max_iteration,1);
% Convergence_curve=zeros(runs,Max_iteration);
% Best_position=zeros(runs,dim);
% %display (['Function:   ', num2str(fn)]);
% for run=1:runs
%     %rng('shuffle');
%     [Best_score,Best_pos,cg_curve]=BSLO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%     Best_score_T(run) = Best_score;
% %     Convergence_curve(run,:)=cg_curve;
%     Best_position(run,:)=Best_pos;
% end
%     Best_score_Best = min(Best_score_T);
%     Best_Score_Mean = mean(Best_score_T);
%     Best_Score_std = std(Best_score_T); 
%     %output
%     format long
%     display([Function_name, ' Best:  ', num2str(Best_score_Best), '     ', 'Mean:  ', num2str(Best_Score_Mean), '     ', 'Std. Deviation:  ', num2str(Best_Score_std)]);
% 
%     x=zeros(1,Max_iteration);
%     for tt=1:1000
%         x(tt)=tt;
%     end
%     % plot(x,AvgConvCurve,'Color','b','LineWidth',3);
%     % title('Average fitness','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     % xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     % set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
% 
%     figure;
%     semilogy(x,cg_curve,'Color','b','LineWidth',3);
%     title('Convergence curve','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
%     t=t+1;
% end
% 
% end
% 
% tt_Max=10;
% runs = 1;
% 
% for fn=24:33
% if fn==24
%     tt_Max=10;
% elseif fn==33
%     tt_Max=10;
% else
%     tt_Max=5;
% end
%         if fn == 23
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 25
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 26 || fn == 27 || fn == 28 || fn == 29 || fn == 30 || fn == 31 || fn == 32
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
% t=1;
% while t<=tt_Max
% initial_flag = 0;
% Function_name=strcat('F',num2str(fn));
% [lb,ub,dim,fobj]=Get_cec20_func(Function_name); 
% Max_iteration= 1000;
%    if dim <= 10
%        MaxFES = 1e5;       
%    elseif dim <= 30 ||dim > 10
%        MaxFES = 2e5;
%    elseif dim <= 50 ||dim > 30
%        MaxFES = 4e5;
%    elseif dim <= 150||dim > 50
%        MaxFES = 8e5;
%    else
%        MaxFES = 1e6;
%    end
% SearchAgents_no=MaxFES/Max_iteration;   
% % Calling algorithm
% Best_score_T = zeros(runs,1);
% AvgConvCurve = zeros(Max_iteration,1);
% Convergence_curve=zeros(runs,Max_iteration);
% Best_position=zeros(runs,dim);
% %display (['Function:   ', num2str(fn)]);
% for run=1:runs
%     %rng('shuffle');
%     [Best_score,Best_pos,cg_curve]=HHO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%     Best_score_T(run) = Best_score;
% %     Convergence_curve(run,:)=cg_curve;
%     Best_position(run,:)=Best_pos;
% end
%     Best_score_Best = min(Best_score_T);
%     Best_Score_Mean = mean(Best_score_T);
%     Best_Score_std = std(Best_score_T); 
%     %output
%     format long
%     display([Function_name, ' Best:  ', num2str(Best_score_Best), '     ', 'Mean:  ', num2str(Best_Score_Mean), '     ', 'Std. Deviation:  ', num2str(Best_Score_std)]);
% 
%     x=zeros(1,Max_iteration);
%     for tt=1:1000
%         x(tt)=tt;
%     end
%     % plot(x,AvgConvCurve,'Color','b','LineWidth',3);
%     % title('Average fitness','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     % xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     % set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
% 
%     figure;
%     semilogy(x,cg_curve,'Color','b','LineWidth',3);
%     title('Convergence curve','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
%     t=t+1;
% end
% 
% end
% 
% tt_Max=10;
% runs = 1;
% 
% for fn=24:33
% if fn==24
%     tt_Max=10;
% elseif fn==33
%     tt_Max=10;
% else
%     tt_Max=5;
% end
%         if fn == 23
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 25
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 26 || fn == 27 || fn == 28 || fn == 29 || fn == 30 || fn == 31 || fn == 32
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
% t=1;
% while t<=tt_Max
% initial_flag = 0;
% Function_name=strcat('F',num2str(fn));
% [lb,ub,dim,fobj]=Get_cec20_func(Function_name); 
% Max_iteration= 1000;
%    if dim <= 10
%        MaxFES = 1e5;       
%    elseif dim <= 30 ||dim > 10
%        MaxFES = 2e5;
%    elseif dim <= 50 ||dim > 30
%        MaxFES = 4e5;
%    elseif dim <= 150||dim > 50
%        MaxFES = 8e5;
%    else
%        MaxFES = 1e6;
%    end
% SearchAgents_no=MaxFES/Max_iteration;   
% % Calling algorithm
% Best_score_T = zeros(runs,1);
% AvgConvCurve = zeros(Max_iteration,1);
% Convergence_curve=zeros(runs,Max_iteration);
% Best_position=zeros(runs,dim);
% %display (['Function:   ', num2str(fn)]);
% for run=1:runs
%     %rng('shuffle');
%     [Best_score,Best_pos,cg_curve]=AOA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%     Best_score_T(run) = Best_score;
% %     Convergence_curve(run,:)=cg_curve;
%     Best_position(run,:)=Best_pos;
% end
%     Best_score_Best = min(Best_score_T);
%     Best_Score_Mean = mean(Best_score_T);
%     Best_Score_std = std(Best_score_T); 
%     %output
%     format long
%     display([Function_name, ' Best:  ', num2str(Best_score_Best), '     ', 'Mean:  ', num2str(Best_Score_Mean), '     ', 'Std. Deviation:  ', num2str(Best_Score_std)]);
% 
%     x=zeros(1,Max_iteration);
%     for tt=1:1000
%         x(tt)=tt;
%     end
%     % plot(x,AvgConvCurve,'Color','b','LineWidth',3);
%     % title('Average fitness','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     % xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     % set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
% 
%     figure;
%     semilogy(x,cg_curve,'Color','b','LineWidth',3);
%     title('Convergence curve','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
%     t=t+1;
% end
% 
% end
% 
% tt_Max=10;
% runs = 1;
% 
% for fn=24:33
% if fn==24
%     tt_Max=10;
% elseif fn==33
%     tt_Max=10;
% else
%     tt_Max=5;
% end
%         if fn == 23
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 25
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 26 || fn == 27 || fn == 28 || fn == 29 || fn == 30 || fn == 31 || fn == 32
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
% t=1;
% while t<=tt_Max
% initial_flag = 0;
% Function_name=strcat('F',num2str(fn));
% [lb,ub,dim,fobj]=Get_cec20_func(Function_name); 
% Max_iteration= 1000;
%    if dim <= 10
%        MaxFES = 1e5;       
%    elseif dim <= 30 ||dim > 10
%        MaxFES = 2e5;
%    elseif dim <= 50 ||dim > 30
%        MaxFES = 4e5;
%    elseif dim <= 150||dim > 50
%        MaxFES = 8e5;
%    else
%        MaxFES = 1e6;
%    end
% SearchAgents_no=MaxFES/Max_iteration;   
% % Calling algorithm
% Best_score_T = zeros(runs,1);
% AvgConvCurve = zeros(Max_iteration,1);
% Convergence_curve=zeros(runs,Max_iteration);
% Best_position=zeros(runs,dim);
% %display (['Function:   ', num2str(fn)]);
% for run=1:runs
%     %rng('shuffle');
%     [Best_score,Best_pos,cg_curve]=RSA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%     Best_score_T(run) = Best_score;
% %     Convergence_curve(run,:)=cg_curve;
%     Best_position(run,:)=Best_pos;
% end
%     Best_score_Best = min(Best_score_T);
%     Best_Score_Mean = mean(Best_score_T);
%     Best_Score_std = std(Best_score_T); 
%     %output
%     format long
%     display([Function_name, ' Best:  ', num2str(Best_score_Best), '     ', 'Mean:  ', num2str(Best_Score_Mean), '     ', 'Std. Deviation:  ', num2str(Best_Score_std)]);
% 
%     x=zeros(1,Max_iteration);
%     for tt=1:1000
%         x(tt)=tt;
%     end
%     % plot(x,AvgConvCurve,'Color','b','LineWidth',3);
%     % title('Average fitness','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     % xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     % set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
% 
%     figure;
%     semilogy(x,cg_curve,'Color','b','LineWidth',3);
%     title('Convergence curve','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
%     t=t+1;
% end
% 
% end
% 
% tt_Max=10;
% runs = 1;
% 
% for fn=24:33
% if fn==24
%     tt_Max=10;
% elseif fn==33
%     tt_Max=10;
% else
%     tt_Max=5;
% end
%         if fn == 23
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 25
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 26 || fn == 27 || fn == 28 || fn == 29 || fn == 30 || fn == 31 || fn == 32
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
% t=1;
% while t<=tt_Max
% initial_flag = 0;
% Function_name=strcat('F',num2str(fn));
% [lb,ub,dim,fobj]=Get_cec20_func(Function_name); 
% Max_iteration= 1000;
%    if dim <= 10
%        MaxFES = 1e5;       
%    elseif dim <= 30 ||dim > 10
%        MaxFES = 2e5;
%    elseif dim <= 50 ||dim > 30
%        MaxFES = 4e5;
%    elseif dim <= 150||dim > 50
%        MaxFES = 8e5;
%    else
%        MaxFES = 1e6;
%    end
% SearchAgents_no=MaxFES/Max_iteration;   
% % Calling algorithm
% Best_score_T = zeros(runs,1);
% AvgConvCurve = zeros(Max_iteration,1);
% Convergence_curve=zeros(runs,Max_iteration);
% Best_position=zeros(runs,dim);
% %display (['Function:   ', num2str(fn)]);
% for run=1:runs
%     %rng('shuffle');
%     [Best_score,Best_pos,cg_curve]=GJO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%     Best_score_T(run) = Best_score;
% %     Convergence_curve(run,:)=cg_curve;
%     Best_position(run,:)=Best_pos;
% end
%     Best_score_Best = min(Best_score_T);
%     Best_Score_Mean = mean(Best_score_T);
%     Best_Score_std = std(Best_score_T); 
%     %output
%     format long
%     display([Function_name, ' Best:  ', num2str(Best_score_Best), '     ', 'Mean:  ', num2str(Best_Score_Mean), '     ', 'Std. Deviation:  ', num2str(Best_Score_std)]);
% 
%     x=zeros(1,Max_iteration);
%     for tt=1:1000
%         x(tt)=tt;
%     end
%     % plot(x,AvgConvCurve,'Color','b','LineWidth',3);
%     % title('Average fitness','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     % xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     % set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
% 
%     figure;
%     semilogy(x,cg_curve,'Color','b','LineWidth',3);
%     title('Convergence curve','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
%     t=t+1;
% end
% 
% end
% 
% tt_Max=10;
% runs = 1;
% 
% for fn=24:33
% if fn==24
%     tt_Max=10;
% elseif fn==33
%     tt_Max=10;
% else
%     tt_Max=5;
% end
%         if fn == 23
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 25
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 26 || fn == 27 || fn == 28 || fn == 29 || fn == 30 || fn == 31 || fn == 32
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
% t=1;
% while t<=tt_Max
% initial_flag = 0;
% Function_name=strcat('F',num2str(fn));
% [lb,ub,dim,fobj]=Get_cec20_func(Function_name); 
% Max_iteration= 1000;
%    if dim <= 10
%        MaxFES = 1e5;       
%    elseif dim <= 30 ||dim > 10
%        MaxFES = 2e5;
%    elseif dim <= 50 ||dim > 30
%        MaxFES = 4e5;
%    elseif dim <= 150||dim > 50
%        MaxFES = 8e5;
%    else
%        MaxFES = 1e6;
%    end
% SearchAgents_no=MaxFES/Max_iteration;   
% % Calling algorithm
% Best_score_T = zeros(runs,1);
% AvgConvCurve = zeros(Max_iteration,1);
% Convergence_curve=zeros(runs,Max_iteration);
% Best_position=zeros(runs,dim);
% %display (['Function:   ', num2str(fn)]);
% for run=1:runs
%     %rng('shuffle');
%     [Best_score,Best_pos,cg_curve]=HOA_v2(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%     Best_score_T(run) = Best_score;
% %     Convergence_curve(run,:)=cg_curve;
%     Best_position(run,:)=Best_pos;
% end
%     Best_score_Best = min(Best_score_T);
%     Best_Score_Mean = mean(Best_score_T);
%     Best_Score_std = std(Best_score_T); 
%     %output
%     format long
%     display([Function_name, ' Best:  ', num2str(Best_score_Best), '     ', 'Mean:  ', num2str(Best_Score_Mean), '     ', 'Std. Deviation:  ', num2str(Best_Score_std)]);
% 
%     x=zeros(1,Max_iteration);
%     for tt=1:1000
%         x(tt)=tt;
%     end
%     % plot(x,AvgConvCurve,'Color','b','LineWidth',3);
%     % title('Average fitness','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     % xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     % set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
% 
%     figure;
%     semilogy(x,cg_curve,'Color','b','LineWidth',3);
%     title('Convergence curve','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
%     t=t+1;
% end
% 
% end
% 
% tt_Max=10;
% runs = 1;
% 
% for fn=24:33
% if fn==24
%     tt_Max=10;
% elseif fn==33
%     tt_Max=10;
% else
%     tt_Max=5;
% end
%         if fn == 23
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 25
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
%         if fn == 26 || fn == 27 || fn == 28 || fn == 29 || fn == 30 || fn == 31 || fn == 32
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
% t=1;
% while t<=tt_Max
% initial_flag = 0;
% Function_name=strcat('F',num2str(fn));
% [lb,ub,dim,fobj]=Get_cec20_func(Function_name); 
% Max_iteration= 1000;
%    if dim <= 10
%        MaxFES = 1e5;       
%    elseif dim <= 30 ||dim > 10
%        MaxFES = 2e5;
%    elseif dim <= 50 ||dim > 30
%        MaxFES = 4e5;
%    elseif dim <= 150||dim > 50
%        MaxFES = 8e5;
%    else
%        MaxFES = 1e6;
%    end
% SearchAgents_no=MaxFES/Max_iteration;   
% % Calling algorithm
% Best_score_T = zeros(runs,1);
% AvgConvCurve = zeros(Max_iteration,1);
% Convergence_curve=zeros(runs,Max_iteration);
% Best_position=zeros(runs,dim);
% %display (['Function:   ', num2str(fn)]);
% for run=1:runs
%     %rng('shuffle');
%     [Best_score,Best_pos,cg_curve]=PO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
%     Best_score_T(run) = Best_score;
% %     Convergence_curve(run,:)=cg_curve;
%     Best_position(run,:)=Best_pos;
% end
%     Best_score_Best = min(Best_score_T);
%     Best_Score_Mean = mean(Best_score_T);
%     Best_Score_std = std(Best_score_T); 
%     %output
%     format long
%     display([Function_name, ' Best:  ', num2str(Best_score_Best), '     ', 'Mean:  ', num2str(Best_Score_Mean), '     ', 'Std. Deviation:  ', num2str(Best_Score_std)]);
% 
%     x=zeros(1,Max_iteration);
%     for tt=1:1000
%         x(tt)=tt;
%     end
%     % plot(x,AvgConvCurve,'Color','b','LineWidth',3);
%     % title('Average fitness','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     % xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     % set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
% 
%     figure;
%     semilogy(x,cg_curve,'Color','b','LineWidth',3);
%     title('Convergence curve','FontSize',21,'FontWeight','bold','FontName','Times New Roman')
%     xlabel('Iteration','FontSize',21,'FontName','Times New Roman');
%     set(gca,'xtick',0:200:1000,'FontName','Times New Roman','FontSize',21);
%     t=t+1;
% end
% 
% end