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
function [Leeches_best_score,Leeches_best_pos,Convergence_curve,sh_x,sh_y]=BSLO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
%% initialize best Leeches
Leeches_best_pos=zeros(1,dim);
Leeches_best_score=inf; 
%search history
sh_x=zeros(1,SearchAgents_no*Max_iter);
sh_y=zeros(1,SearchAgents_no*Max_iter);
%% Initialize the positions of search agents
Leeches_Positions=initialization(SearchAgents_no,dim,ub,lb);
Convergence_curve=zeros(1,Max_iter);
Temp_best_fitness=zeros(1,Max_iter);
fitness=zeros(1,SearchAgents_no);
%% Initialize parameters
t=1;
m=0.8;a=0.97;b=0.001;b2=0.00001;t1=20;t2=20;
% Main loop
while t<=Max_iter
    N1=floor((m+(1-m)*(t/Max_iter)^2)*SearchAgents_no);
    %calculate fitness values
     for i=1:size(Leeches_Positions,1)  
         % boundary checking
         Flag4ub=Leeches_Positions(i,:)>ub;
         Flag4lb=Leeches_Positions(i,:)<lb;
         Leeches_Positions(i,:)=(Leeches_Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
         % Calculate objective function for each search agent
         fitness(i)=fobj(Leeches_Positions(i,:));
           %draw search history
            sh_x(1,(30*(t-1)+i))=Leeches_Positions(i,1);
            sh_y(1,(30*(t-1)+i))=Leeches_Positions(i,2);
         % Update best Leeches
         if fitness(i)<=Leeches_best_score 
             Leeches_best_score=fitness(i); 
             Leeches_best_pos=Leeches_Positions(i,:);
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
        s=8-1*(-(t/Max_iter)^2+1);
    else
        s=8-7*(-(t/Max_iter)^2+1);
    end 
    beta=-0.5*(t/Max_iter)^6+(t/Max_iter)^4+1.5;
    LV=0.5*levy(SearchAgents_no,dim,beta);
    %% Generate random integers
    minValue = 1;  % minimum integer value
    maxValue = floor(SearchAgents_no*(1+t/Max_iter)); % maximum integer value
    k2 = randi([minValue, maxValue], SearchAgents_no, dim);
    k = randi([minValue, dim], SearchAgents_no, dim);
    for i=1:N1
        for j=1:size(Leeches_Positions,2) 
            r1=2*rand()-1; % r1 is a random number in [0,1]
            r2=2*rand()-1;
            r3=2*rand()-1;           
            PD=s*(1-(t/Max_iter))*r1;
            if abs(PD)>=1
                % Exploration of directional leeches
                % b=0.001;
                W1=(1-t/Max_iter)*b*LV(i,j);
                L1=r2*abs(Prey_Position(j)-Leeches_Positions(i,j))*PD*(1-k2(i,j)/SearchAgents_no);
                L2=abs(Prey_Position(j)-Leeches_Positions(i,k(i,j)))*PD*(1-(r2^2)*(k2(i,j)/SearchAgents_no));
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
                % Exploitation of directional leeches
                if t>=0.1*Max_iter
                    b1=b2;
                else
                    b1=b;
                end
                W1=(1-t/Max_iter)*b1*LV(i,j);
                L3=abs(Prey_Position(j)-Leeches_Positions(i,j))*PD*(1-r3*k2(i,j)/SearchAgents_no);
                L4=abs(Prey_Position(j)-Leeches_Positions(i,k(i,j)))*PD*(1-r3*k2(i,j)/SearchAgents_no);
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
        end
    end
    %Search strategy of directionless Leeches
    for i=(N1+1):size(Leeches_Positions,1)
        for j=1:size(Leeches_Positions,2)
            if min(lb)>=0
                LV(i,j)=abs(LV(i,j));
            end
            if rand()>0.5
            Leeches_Positions(i,j)=(t/Max_iter)*LV(i,j)*Leeches_Positions(i,j)*abs(Prey_Position(j)-Leeches_Positions(i,j));  
            else
            Leeches_Positions(i,j)=(t/Max_iter)*LV(i,j)*Prey_Position(j)*abs(Prey_Position(j)-Leeches_Positions(i,j));          
            end
        end
    end
       
    Convergence_curve(t)=Leeches_best_score;   
    t=t+1; 
end
