
%%
clear;
clc
warning off;
resultdir = 'results/randmiss/';
if (~exist('results/randmiss/','file'))
    mkdir('results/randmiss/');
    addpath(genpath('results/randmiss/'));
end
datadir = '../../setting/random remove/';
dataname = {'ORL','orlRnSp','buaaRnSp','caltech7','mfeatRnSp','100Leaves'}; 

datanum = length(dataname);

for datai = 1:6
    datafile = [datadir, cell2mat(dataname(datai))];
    load(datafile); % truth,data,per10-per70,
    fprintf('%s...\n',datafile);

    num_clusters = length(unique(truth));
    num_views    = length(data);
    label = truth;
    r = num_clusters;
    clear truth
    
    for per_in = 1:5  % per incomplete ratio 
        in_ratio = per_in*10;
        percent = per{per_in};%percent = cell(1,10);
        disp(['random miss ', num2str(in_ratio)]);

        for folds = 1:1
            foldspath = [resultdir, char(dataname(datai)),'/',num2str(in_ratio)];
            if (~exist(foldspath,'file'))
                mkdir(foldspath);
                addpath(genpath([foldspath,'/']));
            end 

            index = percent{folds};
            savetxt = [resultdir ,'randmiss_',char(dataname(datai)),'_',num2str(in_ratio),'%','.txt'];
            mes = [char(datetime),'   Folds = ', num2str(folds)];
            dlmwrite(savetxt, mes,'-append','delimiter','\t','newline','pc');

            for vi=1:num_views
                item = index(:,vi);
                W{vi} = diag(item);
                temp = find(item==0);

                X{vi} = data{vi}'; % data:d*n
                X{vi}(temp,:)= 0;
                X{vi} = mapminmax(X{vi},0,1);
                X{vi} = NormalizeFea(X{vi},1);% X{vi}:n*d
                X{vi} = X{vi}';
            end

            %%
            temp = [1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5];
            
            for i = 1:length(temp)
                options.afa = temp(i);
                for j = 1:length(temp)
                    options.beta = temp(j);
                    for aa = 1:length(temp)
                        options.lmd1 = temp(aa);
                        for bb = 1:length(temp)
                            options.lmd2 = temp(bb);
                            disp([options.afa,options.beta, options.lmd1, options.lmd2]);
                            for t = 1:1
                            % tic
                            disp(' ');
                            disp(['t=',num2str(t)]);
                            [U, V, B, obj, S, itr_ac, itr_nmi, itr_fscore] = Graph_DAIMC(X,W,label,r,num_views,options);

                            [Ypred, ~] = SpectralClustering(S,r);
                            
                            Ypred = bestMap(label, Ypred);
                            [ac(t), nmi_value(t), fscore(t)] = CalcMetrics(label, Ypred);
                            fprintf('ac: %0.4f \tnmi:%0.4f\t Fscore: %0.4f\t', ac(t), nmi_value(t), fscore(t));
                            
                            Final_results = [options.afa, options.beta,options.lmd1, options.lmd2,...
                                                mean(ac(i,j,t),3), mean(nmi_value(i,j,t),3), mean(fscore(i,j,t),3)];
                            % toc
                            end
                            matname = [resultdir,char(dataname(datai)),'/',num2str(in_ratio),'/',...
                                                num2str(options.afa),'_',num2str(options.beta),'_',...
                                                num2str(options.lmd1),'_',num2str(options.lmd2),'.mat'];
                            save(matname, 'U', 'V', 'B', 'obj','S',...
                                                    'ac', 'nmi_value', 'fscore','itr_ac', 'itr_nmi', 'itr_fscore');
                            dlmwrite(savetxt, Final_results ,'-append','delimiter','\t','newline','pc');
                        end
                    end
                end
            end
        end
    end
end