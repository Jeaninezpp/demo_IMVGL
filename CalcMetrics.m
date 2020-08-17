function [AC, nmi_value, f_score] = CalcMetrics(label, result)
% Written by Jialu Liu
result = bestMap(label, result);
error_cnt = sum(label ~= result);
AC = length(find(label == result))/length(label);

nmi_value = nmi(label, result);

f_score = compute_f(label, result);