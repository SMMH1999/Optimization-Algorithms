%_________________________________________________________________________________
%  Ameliorated Golden jackal optimization (AGJO)                                                                    
%                                                                                                     
%  Developed in MATLAB R2022a                                                                  
%                                                                                                     
%  programming: Jianfu Bai                                                          
%                                                                                                     
%  e-Mail: Jianfu.Bai@UGent.be, magd.abdelwahab@ugent.be                                                               
%  Soete Laboratory, Department of Electrical Energy, Metals, Mechanical Constructions, and Systems, 
%  Faculty of Engineering and Architecture, Ghent University, Belgium                                                           
%                                                                                                                                                                                                                                                              
%  paper: Jianfu Bai, Magd Abdel Wahab, Ameliorated Golden jackal optimization (AGJO) (2023).  
%____________________________________________________________________________________

function [Male_Jackal_score,Male_Jackal_pos,Convergence_curve]=AGJO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
%% initialize Golden jackal pair
Male_Jackal_pos=zeros(1,dim);
Male_Jackal_score=inf; 
Female_Jackal_pos=zeros(1,dim);  
Female_Jackal_score=inf; 

%% Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
Convergence_curve=zeros(1,Max_iter);
%% Initialize parameters
a=0.1;
b=0.12;
c=0.02;
d=0.8;
t=0;% Loop counter
% Main loop
while t<Max_iter
     for i=1:size(Positions,1)  
         % boundary checking
         Flag4ub=Positions(i,:)>ub;
         Flag4lb=Positions(i,:)<lb;
         Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
         % Calculate objective function for each search agent
         fitness=fobj(Positions(i,:));

         % Update Male Jackal 
         if fitness<Male_Jackal_score 
             Male_Jackal_score=fitness; 
             Male_Jackal_pos=Positions(i,:);
         end  
         if fitness>Male_Jackal_score && fitness<Female_Jackal_score 
             Female_Jackal_score=fitness; 
             Female_Jackal_pos=Positions(i,:);
         end
    end
    E1=1.5*(1-(t/Max_iter));
    %enhanced movement strategy  
    beta=1.5*exp(log(4/3)*t/Max_iter);%beta will be increased from 1.5 to 2
    RL=1.5*levy(SearchAgents_no,dim,beta);
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)  
            r1=rand(); % r1 is a random number in [0,1]
            E0=2*r1-1; 
            E=E1*E0;
            if abs(E)<1
                %% EXPLOITATION
                D_male_jackal=RL(i,j)*abs((Male_Jackal_pos(j)-Positions(i,j))); 
                Male_Positions(i,j)=Male_Jackal_pos(j)+E*D_male_jackal;
                D_female_jackal=RL(i,j)*abs((Female_Jackal_pos(j)-Positions(i,j))); 
                Female_Positions(i,j)=Female_Jackal_pos(j)+E*D_female_jackal;
                
            else
               if t>floor(a*Max_iter)
                %% EXPLORATION
                  D_male_jackal=0.5*RL(i,j)*abs((Male_Jackal_pos(j)-Positions(i,j)));
                  Male_Positions(i,j)=Male_Jackal_pos(j)+E*D_male_jackal;
                  D_female_jackal=0.5*RL(i,j)*abs((Female_Jackal_pos(j)-Positions(i,j)));
                  Female_Positions(i,j)=Female_Jackal_pos(j)+E*D_female_jackal;
               else
                   %use the global search strategy in the beginning
                  if size(ub,2)==1
                   %% EXPLORATION
                    D_male_jackal=0.1*RL(i,j)*abs((ub-Positions(i,j)));
                    Male_Positions(i,j)=Male_Jackal_pos(j)+E*D_male_jackal;
                    D_female_jackal=0.1*RL(i,j)*abs((Positions(i,j)-lb));
                    Female_Positions(i,j)=Female_Jackal_pos(j)+E*D_female_jackal;
                  else
                  %% EXPLORATION
                    D_male_jackal=0.1*RL(i,j)*abs((ub(j)-Positions(i,j)));
                    Male_Positions(i,j)=Male_Jackal_pos(j)+E*D_male_jackal;
                    D_female_jackal=0.1*RL(i,j)*abs((Positions(i,j)-lb(j)));
                    Female_Positions(i,j)=Female_Jackal_pos(j)+E*D_female_jackal; 
                  end
               end
            end             
            %The multi-angle position updating functions for prey
            r2=rand();% r2 is a random number in [0,1]
            r3=rand();% r3 is a random number in [0,1]
            if r2>b
               if t>floor(c*Max_iter)
                  Positions(i,j)=(Male_Positions(i,j)+Female_Positions(i,j))/2;
               else% add a position uncertainty parameter
                  Positions(i,j)=d*(Male_Positions(i,j)+Female_Positions(i,j))/2;
               end
            else %The prey will escape from jackals
               Positions(i,j)=0.5*r3*RL(i,j)*Positions(i,j);
            end                                
        end
    end
    t=t+1;            
    Convergence_curve(t)=Male_Jackal_score;
end



