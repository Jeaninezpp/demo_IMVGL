function [U,V,B,F,S, ac, nmi_value,  fscore] = Graph_DAIMC(X,W,label,r,viewNum,options)
fprintf(sprintf('Initialization...\n'));
[U,V,B] = init(X,W,label,r,viewNum);
% printResult(V, label, r, 1);
S = V*V'; % init S

num_sample = length(label);
eta = 1e-10;
F = 0;
P = 0;
N = 0;
D = cell(viewNum,1);
for i = 1:viewNum
    for k = 1:size(B{i},1)
        D{i}(k,k) = 1/sqrt(norm(B{i}(k,:),2).^2+eta);    
    end
end

time = 0;
f =0;
while 1
    time = time+1;    
    % Update Ui
    for i = 1:viewNum
        tmp1 = options.afa*B{i}*B{i}';
        tmp2 = V'*W{i}*V;
        tmp3 = X{i}*W{i}*V + options.afa*B{i};
        U{i} = lyap(tmp1, tmp2, -tmp3);
    end  
    
    % Update V
    L = diag(sum(S)) - (S+S')/2;
    V = UpdateV(X,W,U,V,L,viewNum,options.lmd1);
    Q = diag(ones(1,size(V,1))*V);
    V = V * inv(Q);
    
    % Update B
    for i = 1:viewNum
        U{i} = U{i}*Q;
        invD = diag(1./diag(0.5*options.beta*D{i}));
        B{i} = (invD - invD * U{i} * inv(U{i}'*invD*U{i} + eye(r)) * U{i}' * invD)*U{i};
       for k = 1:size(B{i},1)
           D{i}(k,k) = 1/sqrt(norm(B{i}(k,:),2).^2+eta);   
       end
    end
    
    % Update S
    for ii = 1:num_sample
        for jj = 1:num_sample
            Q(ii, jj) = (norm(V(ii,:)-V(jj,:)))^2;
        end
    end
    S = - options.lmd1/(2*options.lmd2) * Q';
    
    S = S - diag(diag(S));
    for in = 1:size(S,1)
        idx = 1:size(S,1);
        idx(in) = [];
        S(in,idx) = EProjSimplex_new(S(in,idx));
    end
    
    S(S<0)=0;
    
%     colsum = sum(S,1)+eps;
%     colsum_diag = diag(colsum);
%     S = S * colsum_diag^-1;
%     
    S = (S+S')/2;
    L = diag(sum(S)) - S;
    
    ff = 0;
    for i = 1:viewNum
        tmp1 = (X{i} - U{i}*V')*W{i};
        tmp2 = B{i}'*U{i} - eye(r);
        tmp3 = sum(1./diag(D{i}));
        ff = ff + sum(sum(tmp1.^2)) + options.afa*(sum(sum(tmp2.^2)) + options.beta*tmp3);
    end
    tmp4 = options.lmd1 * trace(V'*L*V);
    tmp5 = options.lmd2 * norm(S,'fro')^2;
    ff = ff + tmp4 + tmp5;
    F(time) = ff;
    
%     indic1 = litekmeans(V, r, 'Replicates', 20);
%     [ac1, nmi_value1] = CalcMetrics(label, indic1);
%     fprintf('V- Finish the iteration %d, ac is %d, nmi is %d, value is %d\n',time,ac1,nmi_value1,ff);

    [indic, ~] = SpectralClustering(S,r);
    [ac(time), nmi_value(time), fscore(time)] = CalcMetrics(label, indic);
    fprintf('Finish the iteration %d, ac is %d, nmi is %d, value is %d\n', time, ac(time), nmi_value(time), ff);

    if abs(ff-f)/f < 1e-4 | abs(ff-f) >  1e100 | time == 50
        break;
    end
    f = ff;
end