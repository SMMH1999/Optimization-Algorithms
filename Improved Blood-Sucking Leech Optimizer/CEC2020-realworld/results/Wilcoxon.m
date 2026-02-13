
clc
clear all
close all % close all figure dialog

maxFunc=16;
% dataset from excel
file_path = '\IBSLO_RC_CEC2020_results.xlsx';
IBSLO = readmatrix(file_path);
file_path = '\BSLO_RC_CEC2020_results.xlsx';
BSLO = readmatrix(file_path);
file_path = '\HHO_RC_CEC2020_results.xlsx';
HHO = readmatrix(file_path);
file_path = '\AOA_RC_CEC2020_results.xlsx';
AOA = readmatrix(file_path);
file_path = '\RSA_RC_CEC2020_results.xlsx';
RSA = readmatrix(file_path);
file_path = '\GJO_RC_CEC2020_results.xlsx';
GJO = readmatrix(file_path);
file_path = '\HOA_RC_CEC2020_results.xlsx';
HOA = readmatrix(file_path);
file_path = '\PO_RC_CEC2020_results.xlsx';
PO = readmatrix(file_path);

p_wilcoxon=zeros(maxFunc,7);

for fn = 1:maxFunc
%         if fn == 2
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
    p_wilcoxon(fn,1)= ranksum(IBSLO(:,fn),BSLO(:,fn));
end
for fn = 1:maxFunc
%         if fn == 2
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
    p_wilcoxon(fn,2)= ranksum(IBSLO(:,fn),HHO(:,fn));
end
for fn = 1:maxFunc
%         if fn == 2
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
    p_wilcoxon(fn,3)= ranksum(IBSLO(:,fn),AOA(:,fn));
end
for fn = 1:maxFunc
%         if fn == 2
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
    p_wilcoxon(fn,4)= ranksum(IBSLO(:,fn),RSA(:,fn));
end
for fn = 1:maxFunc
%         if fn == 2
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
    p_wilcoxon(fn,5)= ranksum(IBSLO(:,fn),GJO(:,fn));
end
for fn = 1:maxFunc
%         if fn == 2
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
    p_wilcoxon(fn,6)= ranksum(IBSLO(:,fn),HOA(:,fn));
end
for fn = 1:maxFunc
%         if fn == 2
%             continue;   %To skip function-2 of CEC-BC-2017 because of its unstable behavior
%         end
    p_wilcoxon(fn,7)= ranksum(IBSLO(:,fn),PO(:,fn));
end
% 写入Excel文件
filename = 'RC_CEC2020_P_Wilcoxon_results.xlsx';  % 定义Excel文件名
% 检查文件是否存在
if exist(filename, 'file')
    % 如果文件存在，则删除文件
    delete(filename);
    disp(['Deleted existing file: ', filename]);
else
    % 如果文件不存在，显示消息
    disp(['File not found, nothing to delete: ', filename]);
end
writematrix(p_wilcoxon, filename);
