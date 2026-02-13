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
function [Leeches_best_score,Leeches_best_pos,Convergence_curve,sh_x,sh_y,Trajectory,AvgConvCurve]=IBSLO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
%% initialize best Leeches
Leeches_best_pos=zeros(1,dim);
Leeches_best_score=inf; 
%search history
sh_x=zeros(1,SearchAgents_no*Max_iter);
sh_y=zeros(1,SearchAgents_no*Max_iter);
Trajectory=zeros(1,Max_iter);
AvgConvCurve=zeros(1,Max_iter);
%% Initialize the positions of search agents
Leeches_Positions=initialization(SearchAgents_no,dim,ub,lb);
Convergence_curve=zeros(1,Max_iter);
Temp_best_fitness=zeros(1,Max_iter);
fitness=zeros(1,SearchAgents_no);
%% Initialize parameters
t=1;
m=1;a=0.97;b=0.001;b2=0.00001;t1=20;t2=20;
CM=chaos(2,Max_iter,1);
CM2=chaos(9,Max_iter,1);
% Main loop

while t<=Max_iter
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

    Trajectory(t)=Leeches_Positions(1,1);
    AvgConvCurve(t)=mean(fitness);
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
    LV=0.05*levy(SearchAgents_no,dim,beta);%0.06
    %% Generate random integers
    minValue = 1;  % minimum integer value
    maxValue = floor(SearchAgents_no*(1+t/Max_iter)); % maximum integer value
    k2 = randi([minValue, maxValue], SearchAgents_no, dim);
    k = randi([minValue, dim], SearchAgents_no, dim);

    for i=1:size(Leeches_Positions,1) 
        k3=randi(dim);
        for j=1:size(Leeches_Positions,2) 
            r1=2*rand()-1; % r1 is a random number in [0,1]
            r2=2*rand()-1;
            r3=2*rand()-1;
            PD=s*(1-(t/Max_iter))*r1;
            PD1=PD*0.8*CM(t);

            if abs(PD1)>=1
                % Exploration of directional leeches
                PD=PD1;
                if i<0.5*size(Leeches_Positions,1)
                     b=1;
                else
                    b=0.001;
                end
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
                if rand()<0.2 && t<=0.8*Max_iter
                    PD=PD1;
                else
                    PD=PD*(2-3*t/Max_iter)*CM2(t);
                end
                % Exploitation of directional leeches
                % if t>=0.1*Max_iter && rand()<0.7
                %     b1=b2;
                % else
                %     if rand()<0.8
                %         b1=1;
                %     else
                %         b1=0.001;
                %     end
                % end
                if i<=0.6*size(Leeches_Positions,1) 
                     b1=b2;
                else
                if t>=0.1*Max_iter && rand()<0.7
                    b1=b2;
                else
                    if rand()<0.8
                        b1=1;
                    else
                        b1=0.001;
                    end
                end
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
            r4=rand();
            if i>0.6*size(Leeches_Positions,1) && t<0.4*Max_iter 
                if size(ub,2)<=1
                    Leeches_Positions(i,j)=(1-t/Max_iter)*(ub+lb)-Leeches_Positions(i,k3);
                else
                    Leeches_Positions(i,j)=(1-t/Max_iter)*(ub(j)+lb(j))-Leeches_Positions(i,k3);
                end
            elseif r4>0.2 && t<0.1*Max_iter
                    Leeches_Positions(i,j)=(t/Max_iter)*LV(i,j)*Leeches_Positions(i,j)*abs(Prey_Position(j)-Leeches_Positions(i,j));           
            elseif i>=(m+0.2*(t/Max_iter)^2)*size(Leeches_Positions,1) && t>=0.1*Max_iter 
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
    end

    Convergence_curve(t)=Leeches_best_score;   
    t=t+1; 
end
