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
% This function creates the first random population

function X=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal
if Boundary_no==1
    X=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end