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

clc
clear all
close all % close all figure dialog
tic
disp('enter ANN_IBSLO');
% Define the neural network architecture
input_nodes = 5;
hidden_nodes = 14;
output_nodes = 1;

% training dataset from excel
file_path = '\NFES_Data1_2024.xlsx';
M = readmatrix(file_path);
X_train=M(:,1:5);
y_train=M(:,6);

% % 归一化训练特征
% X_min = min(X_train);
% X_max = max(X_train);
% X_normalized = (X_train - X_min) ./ (X_max - X_min);
% 
% % 保存目标变量的最小值和最大值，并归一化目标变量
% y_min = min(y_train);
% y_max = max(y_train);
% y_normalized = (y_train - y_min) / (y_max - y_min);
% 标准化 X_train
X_mean = mean(X_train);
X_std = std(X_train);
X_train_normalized = (X_train - X_mean) ./ X_std;

% 如果 y_train 也需要归一化（取决于模型是否对目标变量的范围敏感）
y_mean = mean(y_train);
y_std = std(y_train);
y_train_normalized = (y_train - y_mean) ./ y_std;
ii=0;
while ii<=0
% Initialize parameters
run_max=1;% running time
run=1;
Max_iterations = 1000;
population_size = 150;
search_space_min = -1;
search_space_max = 1;
lb=search_space_min;ub=search_space_max;
t=1;m=1;a=0.97;b=0.001;b2=0.00001;t1=20;t2=20;
% initialize the golden jackals
Leeches_best_pos=zeros(1,(input_nodes*hidden_nodes+2*hidden_nodes+1)); 
Temp_best_fitness=zeros(1,Max_iterations);
Convergence_curve=zeros(run_max,Max_iterations);
All_MSE=zeros(1,run_max);
All_Pos=zeros(run_max,(input_nodes*hidden_nodes+2*hidden_nodes+1));
dim=(input_nodes*hidden_nodes+2*hidden_nodes+1);
CM=chaos(2,Max_iterations,1);
CM2=chaos(9,Max_iterations,1);
while run<=run_max
% Step 1: Initialize the preys' population with random weights and biases
Leeches_Positions=search_space_min + (search_space_max - search_space_min) * rand(population_size, (input_nodes*hidden_nodes+2*hidden_nodes+1));%one output
% Main loop
for t = 1:Max_iterations
    N1=floor((m+(1-m)*(t/Max_iterations)^2)*population_size);
    % Step 2: Evaluate fitness (MSE) for each prey in the population
    fitness = zeros(1, population_size);
    for i = 1:population_size
        % boundary checking
        Flag4search_space_max=Leeches_Positions(i,:)>search_space_max;
        Flag4search_space_min=Leeches_Positions(i,:)<search_space_min;
        Leeches_Positions(i,:)=(Leeches_Positions(i,:).*(~(Flag4search_space_max+Flag4search_space_min)))+search_space_max.*Flag4search_space_max+search_space_min.*Flag4search_space_min;    
        weights_input_hidden=Leeches_Positions(i,1:(input_nodes*hidden_nodes));
        weights_hidden_output=Leeches_Positions(i,(input_nodes*hidden_nodes+1):(input_nodes*hidden_nodes+hidden_nodes));
        biases_hidden=Leeches_Positions(i,(input_nodes*hidden_nodes+hidden_nodes+1):(input_nodes*hidden_nodes+2*hidden_nodes));
        biases_output=Leeches_Positions(i,(input_nodes*hidden_nodes+2*hidden_nodes+1));
  
    mse_accumulator = 0;
    % LOOCV loop
    for j = 1:size(X_train_normalized, 1)
        % Create training and test sets by leaving out the j-th sample
        X_train_loocv = X_train_normalized([1:j-1, j+1:end], :);
        y_train_loocv = y_train_normalized([1:j-1, j+1:end]);

        % Prediction for the j-th sample
        X_test_loocv = X_train_normalized(j, :);
        y_test_loocv = y_train_normalized(j);

        % Prediction using the current set of weights and biases
        y_pred_normalized_loocv = neural_network_predict(X_test_loocv, weights_input_hidden, weights_hidden_output, biases_hidden, biases_output);

        % Calculate the error for the j-th sample
        error_loocv = y_test_loocv - y_pred_normalized_loocv;
        mse_accumulator = mse_accumulator + (error_loocv^2);
    end

    % Average the MSE over all samples
    fitness(i) = mse_accumulator / size(X_train_normalized, 1);
    end    
    % Step 3: Update the male and femal jackals (best solution)
    [fitness_sorted, sort_idx]=sort(fitness);
    sorted_population=Leeches_Positions(sort_idx,:);
    min_fitness=fitness_sorted(1);
    leader=sorted_population(1,:);
    if t==1
        Leeches_best_score=min_fitness;
        Leeches_best_pos=leader;
    else
        if min_fitness<Leeches_best_score
            Leeches_best_score=min_fitness;
            Leeches_best_pos=leader;
        end
    end

    Prey_Position=Leeches_best_pos;
     % Re-tracking strategy
     Temp_best_fitness(t)=Leeches_best_score;
     if t>t1
         if Temp_best_fitness(t)==Temp_best_fitness(t-t2)
             for i=1:size(Leeches_Positions,1) 
                 if fitness(i)==Leeches_best_score 
                    Leeches_Positions(i,:)=rand(1,dim).*(ub-lb)+lb;
                 end
             end
         end
     end    

    if rand()<0.5
        s=8-1*(-(t/Max_iterations)^2+1);
    else
        s=8-7*(-(t/Max_iterations)^2+1);
    end 
    beta=-0.5*(t/Max_iterations)^6+(t/Max_iterations)^4+1.5;
    LV=0.05*levy(population_size,dim,beta);%0.06
    %% Generate random integers
    minValue = 1;  % minimum integer value
    maxValue = floor(population_size*(1+t/Max_iterations)); % maximum integer value
    k2 = randi([minValue, maxValue], population_size, dim);
    k = randi([minValue, dim], population_size, dim);
    k3=randi(dim);
    % mo= mean(Leeches_Positions);
    for i=1:N1
        k4=randi(dim);
        for j=1:size(Leeches_Positions,2) 
            r1=2*rand()-1; % r1 is a random number in [0,1]
            r2=2*rand()-1;
            r3=2*rand()-1;           

            PD=s*(1-(t/Max_iterations))*r1;
            PD1=PD*0.8*CM(t);
            if abs(PD1)>=1
                % Exploration of directional leeches
                % b=0.001;
                PD=PD*0.8*CM(t);
                if i<0.5*size(Leeches_Positions,1)
                     b=1;
                else
                    b=0.001;
                end             
                W1=(1-t/Max_iterations)*b*LV(i,j);
                L1=r2*abs(Prey_Position(j)-Leeches_Positions(i,j))*PD*(1-k2(i,j)/population_size);
                L2=abs(Prey_Position(j)-Leeches_Positions(i,k(i,j)))*PD*(1-(r2^2)*(k2(i,j)/population_size));
                if rand()<a
                if abs(Prey_Position(j))>abs(Leeches_Positions(i,j))
                Leeches_Positions(i,j)=Leeches_Positions(i,j)+W1*Leeches_Positions(i,j)-L1;
                else
                Leeches_Positions(i,j)=Leeches_Positions(i,j)+W1*Leeches_Positions(i,j)+L1;
                end
                else
                if abs(Prey_Position(j))>abs(Leeches_Positions(i,j))
                Leeches_Positions(i,j)=Leeches_Positions(i,j)+W1*Leeches_Positions(i,k(i,j))-L2;
                else
                Leeches_Positions(i,j)=Leeches_Positions(i,j)+W1*Leeches_Positions(i,k(i,j))+L2;
                end
                end  
            else
                if rand()<0.2 && t<=0.8*Max_iterations
                PD=PD*0.8*CM(t);
                else
                    % PD=PD*(1.5-1.5*t/Max_iterations)*CM2(t);
                    PD=PD*(2-3*t/Max_iterations)*CM2(t);                    
                end
                % Exploitation of directional leeches
                % if t>=0.1*Max_iterations && rand()<0.7
                %     b1=b2;
                % else
                %     if rand()<0.8
                %     b1=1;
                %     else
                %         b1=0.001;
                %     end
                % end
                if i<=0.6*size(Leeches_Positions,1) 
                     b1=b2;
                else
                if t>=0.1*Max_iterations && rand()<0.7
                    b1=b2;
                else
                    if rand()<0.8
                    b1=1;
                    else
                        b1=0.001;
                    end
                end
                end
                W1=(1-t/Max_iterations)*b1*LV(i,j);
                L3=abs(Prey_Position(j)-Leeches_Positions(i,j))*PD*(1-r3*k2(i,j)/population_size);
                L4=abs(Prey_Position(j)-Leeches_Positions(i,k(i,j)))*PD*(1-r3*k2(i,j)/population_size);
                if rand()<a
                if abs(Prey_Position(j))>abs(Leeches_Positions(i,j))
                Leeches_Positions(i,j)=Prey_Position(j)+W1*Prey_Position(j)-L3;
                else
                Leeches_Positions(i,j)=Prey_Position(j)+W1*Prey_Position(j)+L3;
                end
                else
                 if abs(Prey_Position(j))>abs(Leeches_Positions(i,j))
                Leeches_Positions(i,j)=Prey_Position(j)+W1*Prey_Position(j)-L4;
                else
                Leeches_Positions(i,j)=Prey_Position(j)+W1*Prey_Position(j)+L4;
                end                   
                end
            end
            r4=rand();
            % rv=0.7-0.7*t/Max_iterations;%这个设置结果还行
            if i>0.6*size(Leeches_Positions,1) && t<0.4*Max_iterations %i>0.5*size(Leeches_Positions,1) && t<rv*Max_iter,rv=0.1则F15最优
                if size(ub,2)<=1
                    Leeches_Positions(i,j)=(ub+lb)/2-Leeches_Positions(i,k4);%重点，对角线
                else
                    Leeches_Positions(i,j)=(ub(j)+lb(j))/2-Leeches_Positions(i,k4);
                end            
            elseif r4>0.2 && t<0.1*Max_iterations 
                % if min(lb)>=0
                %     LV(i,j)=abs(LV(i,j));
                % end
                    Leeches_Positions(i,j)=(t/Max_iterations)*LV(i,j)*Leeches_Positions(i,j)*abs(Prey_Position(j)-Leeches_Positions(i,j));  
            elseif i>=(0.8+0.2*(t/Max_iterations)^2)*size(Leeches_Positions,1) && t>=0.1*Max_iterations %这样高纬度下F4也表现最优
                if min(lb)>=0
                    LV(i,j)=abs(LV(i,j));
                end
                if rand()>0.5
                    Leeches_Positions(i,j)=(t/Max_iterations)*LV(i,j)*Leeches_Positions(i,j)*abs(Prey_Position(j)-Leeches_Positions(i,j));  
                else
                    Leeches_Positions(i,j)=(t/Max_iterations)*LV(i,j)*Prey_Position(j)*abs(Prey_Position(j)-Leeches_Positions(i,j));          
                end                 
            end
        end
    end
       
    Convergence_curve(run,t)=Leeches_best_score;   
end
All_MSE(run)=Leeches_best_score;
All_Pos(run,:)=Leeches_best_pos;
run=run+1;
fprintf('MSE in each run: %.6f\n', Leeches_best_score);
end
Ave_MSE=mean(All_MSE);
[Min_MSE, index]=min(All_MSE);
Best_pos=All_Pos(index,:);
Convergence_curve=Convergence_curve(index,:);
Max_MSE=max(All_MSE);
Std_MSE=std(All_MSE);
fprintf('Last_MSE: %.6f\n', Leeches_best_score);
fprintf('Ave_MSE: %.6f\n', Ave_MSE);
fprintf('Min_MSE: %.6f\n', Min_MSE);
fprintf('Max_MSE: %.6f\n', Max_MSE);
fprintf('Std_MSE: %.6f\n', Std_MSE);
% weights and biases with minimum mse is obtained
trained_weights_input_hidden=Best_pos(1:(input_nodes*hidden_nodes));
trained_weights_hidden_output=Best_pos((input_nodes*hidden_nodes+1):(input_nodes*hidden_nodes+hidden_nodes));
trained_biases_hidden=Best_pos((input_nodes*hidden_nodes+hidden_nodes+1):(input_nodes*hidden_nodes+2*hidden_nodes));
trained_biases_output=Best_pos((input_nodes*hidden_nodes+2*hidden_nodes+1));
%draw convergence curves
x=zeros(1,Max_iterations);
for tt=1:Max_iterations
    x(tt)=tt;
end
plot(x,Convergence_curve,'Color','b','LineWidth',3);
title('Convergence curve','FontSize',18,'FontWeight','bold','FontName','Times New Roman')
xlabel('Iteration','FontSize',18,'FontName','Times New Roman');
set(gca,'xtick',0:100:Max_iterations,'FontName','Times New Roman','FontSize',18);
figure;

% y_normalized = rescale(y, 0, 1);
predicted_output2 = neural_network_predict(X_train_normalized, trained_weights_input_hidden, trained_weights_hidden_output, trained_biases_hidden, trained_biases_output);
predicted_diameter = predicted_output2 * y_std + y_mean;
p=plotregression(y_train,predicted_diameter'); 

disp(ii);
ii=ii+1;
toc

end

