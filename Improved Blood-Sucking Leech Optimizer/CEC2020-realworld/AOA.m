        
function [Best_FF,Best_P,Conv_curve]=AOA(N,M_Iter,LB,UB,Dim,F_obj)
% display('AOA Working');
%Two variables to keep the positions and the fitness value of the best-obtained solution

Best_P=zeros(1,Dim);
Best_FF=inf;
Conv_curve=zeros(1,M_Iter);

%Initialize the positions of solution，20个个体代表二十行，每个个体有十列，代表十个维度
X=initialization(N,Dim,UB,LB);
Xnew=X;
Ffun=zeros(1,size(X,1));% (fitness values)，有二十个个体，所以是一行20列的矩阵，用zeros是为了提高计算效率，提前分配存储空间
Ffun_new=zeros(1,size(Xnew,1));% (fitness values)

MOP_Max=1;
MOP_Min=0.2;
C_Iter=2;
Alpha=5;
Mu=0.499;


for i=1:size(X,1)
    Ffun(1,i)=F_obj(X(i,:));  %Calculate the fitness values of solutions，目标函数算出来的标量
    if Ffun(1,i)<Best_FF    %适应度函数设为目标函数，那么适应度值当然越小越好。如果适应度设为目标函数的倒数，那么适应度越大越好。
        Best_FF=Ffun(1,i);
        Best_P=X(i,:);
    end
end
Conv_curve(1)=Best_FF;    
    

while C_Iter<M_Iter+1  %Main loop
    MOP=1-((C_Iter)^(1/Alpha)/(M_Iter)^(1/Alpha));   % Probability Ratio 
    MOA=MOP_Min+C_Iter*((MOP_Max-MOP_Min)/M_Iter); %Accelerated function
   
    %Update the Position of solutions
    for i=1:size(X,1)   % if each of the UB and LB has a just value 
        for j=1:size(X,2)
           r1=rand();
            if (size(LB,2)==1)
                if r1<MOA  %In literature， in in the Fig. 1 for Exploration phase r1<MOA but in Algorithm and flow chart r1>MOA, and there, r1<MOA why?
                    r2=rand();
                    if r2>0.5  % why is r2>0.5 instead of r2<0.5
                        Xnew(i,j)=Best_P(1,j)/(MOP+eps)*((UB-LB)*Mu+LB);
                    else
                        Xnew(i,j)=Best_P(1,j)*MOP*((UB-LB)*Mu+LB);
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,j)=Best_P(1,j)-MOP*((UB-LB)*Mu+LB);
                    else
                        Xnew(i,j)=Best_P(1,j)+MOP*((UB-LB)*Mu+LB);
                    end
                end               
            end
            
           
            if (size(LB,2)~=1)   % if each of the UB and LB has more than one value 
                r1=rand();
                if r1<MOA
                    r2=rand();
                    if r2>0.5
                        Xnew(i,j)=Best_P(1,j)/(MOP+eps)*((UB(j)-LB(j))*Mu+LB(j));
                    else
                        Xnew(i,j)=Best_P(1,j)*MOP*((UB(j)-LB(j))*Mu+LB(j));
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,j)=Best_P(1,j)-MOP*((UB(j)-LB(j))*Mu+LB(j));
                    else
                        Xnew(i,j)=Best_P(1,j)+MOP*((UB(j)-LB(j))*Mu+LB(j));
                    end
                end               
            end
            
        end
        
        Flag_UB=Xnew(i,:)>UB; % check if they exceed (up) the boundaries
        Flag_LB=Xnew(i,:)<LB; % check if they exceed (down) the boundaries
        Xnew(i,:)=(Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB;  %??????
 
        Ffun_new(1,i)=F_obj(Xnew(i,:));  % calculate Fitness function 
        if Ffun_new(1,i)<Ffun(1,i)
            X(i,:)=Xnew(i,:);
            Ffun(1,i)=Ffun_new(1,i);
        end
        if Ffun(1,i)<Best_FF
        Best_FF=Ffun(1,i);
        Best_P=X(i,:);
        end
       
    end
    

    %Update the convergence curve
    Conv_curve(C_Iter)=Best_FF;
    
    %Print the best solution details after every 50 iterations
%     if mod(C_Iter,50)==0
%         display(['At iteration ', num2str(C_Iter), ' the best solution fitness is ', num2str(Best_FF)]);
%     end
     
    C_Iter=C_Iter+1;  % incremental iteration
   
end



