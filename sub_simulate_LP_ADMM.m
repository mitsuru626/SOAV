%% read problem instances from MAT-file and solve them
clearvars -except LoadFileName SaveFileName SaveFileNameList SuffixList gamma n_trial n_trial_conventional ConsiderQuadraticCost
clc;

load(LoadFileName);
n_point = length(nwTable);

% arrays of computation time
% measured by MATLAB tic-toc
mfile_time_bisec_fixed = zeros(n_point,1);
mfile_time_bisec_adaptive = zeros(n_point,1);
mfile_time_interp = zeros(n_point,1);
mfile_time_conventional = zeros(n_point,1);
mfile_time_LP_interior = zeros(n_point,1);
mfile_time_LP_psimplex = zeros(n_point,1);
mfile_time_LP_dsimplex = zeros(n_point,1);
% measured in MEX code
solver_time_bisec_fixed = zeros(n_point,1);
solver_time_bisec_adaptive = zeros(n_point,1);
solver_time_interp = zeros(n_point,1);
solver_time_conventional = zeros(n_point,1);
solver_time_LP_interior = zeros(n_point,1);
solver_time_LP_psimplex = zeros(n_point,1);
solver_time_LP_dsimplex = zeros(n_point,1);

% arrays of numbers of iterations
iter_fixed = zeros(n_point,1);
iter_adaptive = zeros(n_point,1);
iter_interp = zeros(n_point,1);
iter_conventional = zeros(n_point,1);

for i_trial=1:n_trial
    for i_point = 1:n_point
        load(LoadFileName);
        n_w = nwTable(i_point);
        N = NTable(i_point);
        ind = find(n_w==nwTable);
        fprintf('N = %d, n_w = %d: ',N,n_w);
        Aeq = AeqTable{ind,1}; beq = beqTable{ind,1}; n_b = length(beq);
        lb = lbTable{ind,1}; ub = ubTable{ind,1};
        u = uTable{ind,1};
        w = wTable{ind,1};
        y0 = zeros(N,1); z0 = zeros(N,1);
        Q = QTable{ind,1};
        c = cTable{ind,1};
        clearvars AeqTable beqTable lbTable ubTable uTable QTable cTable;
        
        % LP primal simplex
        options = struct(); options.Method = 0; [Uopt_LP_psimplex,fval_LP_psimplex,info_LP_psimplex] = ...
            soav_LP_QP(u,w,lb,ub,Aeq,beq,Q,c,options);
        mfile_time_LP_psimplex(i_point,1) = info_LP_psimplex.MfileTime;
        solver_time_LP_psimplex(i_point,1) = info_LP_psimplex.GurobiTime;
        % LP dual simplex
        options = struct(); options.Method = 1; [Uopt_LP_dsimplex,fval_LP_dsimplex,info_LP_dsimplex] = ...
            soav_LP_QP(u,w,lb,ub,Aeq,beq,Q,c,options);
        mfile_time_LP_dsimplex(i_point,1) = info_LP_dsimplex.MfileTime;
        solver_time_LP_dsimplex(i_point,1) = info_LP_dsimplex.GurobiTime;
        % LP interior
        options = struct(); options.Method = 2; [Uopt_LP_interior,fval_LP_interior,info_LP_interior] = ...
            soav_LP_QP(u,w,lb,ub,Aeq,beq,Q,c,options);
        mfile_time_LP_interior(i_point,1) = info_LP_interior.MfileTime;
        solver_time_LP_interior(i_point,1) = info_LP_interior.GurobiTime;
        % Tolerance Parameters
        ConstraintTol = min(info_LP_interior.ConstraintViolationList);
        fvalTol = min(abs(info_LP_interior.fvalList-fval_LP_psimplex));
        fvalTrue = fval_LP_psimplex;
        fprintf('fvalTol %g, ConstraintTol = %g \n',fvalTol,ConstraintTol);
        MexOptions = [fvalTrue,fvalTol,ConstraintTol];
        
        
        % Proposed; Fixed
        tic;
        [Uopt_bisec_fixed,fval_bisec_fixed,exitflag_bisec_fixed,output_bisec_fixed] = ...
            soav_bisec_fixed(u,w,lb,ub,Aeq,beq,Q,c,gamma,y0,z0,MexOptions);
        if exitflag_bisec_fixed == 1
            mfile_time_bisec_fixed(i_point,1) = toc;
            solver_time_bisec_fixed(i_point,1) =  output_bisec_fixed(1);
            iter_fixed(i_point,1) = output_bisec_fixed(2);
        end
        % Proposed: ADAPTIVE
        tic;
        [Uopt_bisec_adaptive,fval_bisec_adaptive,exitflag_bisec_adaptive,output_bisec_adaptive] = ...
            soav_bisec_adaptive(u,w,lb,ub,Aeq,beq,Q,c,gamma,y0,z0,MexOptions);
        if exitflag_bisec_adaptive == 1
            mfile_time_bisec_adaptive(i_point,1) = toc;
            solver_time_bisec_adaptive(i_point,1) =  output_bisec_adaptive(1);
            iter_adaptive(i_point,1) = output_bisec_adaptive(2);
        end
        
        % Proposed: interp
        tic;
        [Uopt_interp,fval_interp,exitflag_interp,output_interp] = ...
            soav_interp(u,w,lb,ub,Aeq,beq,Q,c,gamma,y0,z0,MexOptions);
        if exitflag_interp == 1
            mfile_time_interp(i_point,1) = toc;
            solver_time_interp(i_point,1) =  output_interp(1);
            iter_interp(i_point,1) = output_interp(2);
        end
        
        % conventional ADMM
        if i_trial <= n_trial_conventional
            tic;
            y0 = zeros((n_w+1)*N + n_b,1); z0 = zeros((n_w+1)*N + n_b,1);
            [Uopt_conventional,fval_conventional,exitflag_conventional,output_conventional]=...
                soav_conventional(u,w,lb,ub,Aeq,beq,gamma,y0,z0,MexOptions);
            if exitflag_conventional == 1
                mfile_time_conventional(i_point,1) = toc;
                solver_time_conventional(i_point,1) =  output_conventional(1);
                iter_conventional(i_point,1) = output_conventional(2);
            else
                mfile_time_conventional(i_point,1) =  -1;
                solver_time_conventional(i_point,1) =  -1;
                iter_conventional(i_point,1) = -1;   
            end
        end
        
        % display
        if i_trial <= n_trial_conventional
        fprintf('fval: Proposed (fixed, adapt, interp) = (%g, %g, %g), LP (interior, simplex) = (%g, %g), Conventional = %g \n',...
            fval_bisec_fixed,fval_bisec_adaptive,fval_interp,fval_LP_interior,...
            fval_LP_psimplex,fval_conventional);
        fprintf('t_iter: Proposed (%f/%f  %f/%f  %f/%f) LP (%f/%f) Conventional (%f/%f) \n',...
            mfile_time_bisec_fixed(i_point,1),solver_time_bisec_fixed(i_point,1),...
            mfile_time_bisec_adaptive(i_point,1),solver_time_bisec_adaptive(i_point,1),...
            mfile_time_interp(i_point,1),solver_time_interp(i_point,1),...
            mfile_time_LP_interior(i_point,1),solver_time_LP_interior(i_point,1),....
            mfile_time_conventional(i_point,1),solver_time_conventional(i_point,1));
        fprintf('n_iter: Proposed (%d, %d, %d) conventional (%d) \n',...
            iter_fixed(i_point,1),iter_adaptive(i_point,1),iter_interp(i_point,1),iter_conventional(i_point,1));
        else
        fprintf('fval: Proposed (fixed, adapt, interp) = (%g, %g, %g), LP (interior, simplex) = (%g, %g) \n',...
            fval_bisec_fixed,fval_bisec_adaptive,fval_interp,fval_LP_interior,...
            fval_LP_psimplex);
        fprintf('t_iter: Proposed (%f/%f  %f/%f  %f/%f) LP (%f/%f) \n',...
            mfile_time_bisec_fixed(i_point,1),solver_time_bisec_fixed(i_point,1),...
            mfile_time_bisec_adaptive(i_point,1),solver_time_bisec_adaptive(i_point,1),...
            mfile_time_interp(i_point,1),solver_time_interp(i_point,1),...
            mfile_time_LP_interior(i_point,1),solver_time_LP_interior(i_point,1));
        fprintf('n_iter: Proposed (%d, %d, %d) \n',...
            iter_fixed(i_point,1),iter_adaptive(i_point,1),iter_interp(i_point,1));
        end
    end
    % save data 
    if i_trial <= n_trial_conventional
        save(['Result_N_',num2str(min(NTable)),'_',num2str(max(NTable)),...
            '_nw_',num2str(min(nwTable)),'_',num2str(max(nwTable)),...
            '_trial_',num2str(i_trial)],...
            '-regexp','^mfile_time','^solver_time','^iter','nwTable','NTable');
        if i_trial == n_trial_conventional
            % experiment of conventional ADMM has been finished and clear
            % the data 
            clear mfile_time_conventional solver_time_conventional iter_conventional
        end
    else
        save(['Result_N_',num2str(min(NTable)),'_',num2str(max(NTable)),...
            '_nw_',num2str(min(nwTable)),'_',num2str(max(nwTable)),...
            '_trial_',num2str(i_trial)],...
            '-regexp','^mfile_time','^solver_time','^iter','nwTable','NTable');
    end
end

%% save data
clearvars -except nwTable NTable SaveFileName SaveFileNameList SuffixList gamma n_trial n_trial_conventional ConsiderQuadraticCost;
regexp_file = ['Result_N_',num2str(min(NTable)),'_',num2str(max(NTable)),'_nw_',num2str(min(nwTable)),'_',num2str(max(nwTable)),'_trial_*'];
suffix_file = ['_N_',num2str(min(NTable)),'_',num2str(max(NTable)),'_nw_',num2str(min(nwTable)),'_',num2str(max(nwTable))];
sub_MATaverage([regexp_file(1:end-1),'1.mat'],'nwTable',['nwTable',suffix_file]); % read from first file 
sub_MATaverage([regexp_file(1:end-1),'1.mat'],'NTable',['NTable',suffix_file]);
sub_MATaverage(regexp_file,'solver_time_LP_interior',['solver_time_LP_interior',suffix_file]);
sub_MATaverage(regexp_file,'solver_time_LP_psimplex',['solver_time_LP_psimplex',suffix_file]);
sub_MATaverage(regexp_file,'solver_time_LP_dsimplex',['solver_time_LP_dsimplex',suffix_file]);
sub_MATaverage(regexp_file,'solver_time_conventional',['solver_time_conventional',suffix_file]);
sub_MATaverage(regexp_file,'solver_time_bisec_fixed',['solver_time_bisec_fixed',suffix_file]);
sub_MATaverage(regexp_file,'solver_time_bisec_adaptive',['solver_time_bisec_adaptive',suffix_file]);
sub_MATaverage(regexp_file,'solver_time_interp',['solver_time_interp',suffix_file]);
sub_MATaverage(regexp_file,'mfile_time_LP_interior',['mfile_time_LP_interior',suffix_file]);
sub_MATaverage(regexp_file,'mfile_time_LP_psimplex',['mfile_time_LP_psimplex',suffix_file]);
sub_MATaverage(regexp_file,'mfile_time_LP_dsimplex',['mfile_time_LP_dsimplex',suffix_file]);
sub_MATaverage(regexp_file,'mfile_time_conventional',['mfile_time_conventional',suffix_file]);
sub_MATaverage(regexp_file,'mfile_time_bisec_fixed',['mfile_time_bisec_fixed',suffix_file]);
sub_MATaverage(regexp_file,'mfile_time_bisec_adaptive',['mfile_time_bisec_adaptive',suffix_file]);
sub_MATaverage(regexp_file,'mfile_time_interp',['mfile_time_interp',suffix_file]);
sub_MATaverage(regexp_file,'iter_fixed',['iter_fixed',suffix_file]);
sub_MATaverage(regexp_file,'iter_adaptive',['iter_adaptive',suffix_file]);
sub_MATaverage(regexp_file,'iter_interp',['iter_interp',suffix_file]);
sub_MATaverage(regexp_file,'iter_conventional',['iter_conventional',suffix_file]);
save(SaveFileName,'-regexp','^mfile_time','^solver_time','^iter','nwTable','NTable');
delete(regexp_file);