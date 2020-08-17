function [ac, nmi_value, fscore] = printResult(X, label, K, kmeansFlag)

if kmeansFlag == 1
    indic = litekmeans(X, K, 'Replicates', 20);
else
    [~, indic] = max(X, [] ,2);
end
indic = bestMap(label, indic);
[ac, nmi_value, fscore] = CalcMetrics(label, indic);
fprintf('ac: %0.4f \tnmi:%0.4f\t Fscore: %0.4f\t', ac, nmi_value, fscore);