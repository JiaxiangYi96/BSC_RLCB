% this code is used to estimate the failure probability using the proposed RLCB+BSC method
% code by Jiaxiang Yi
% Email: jiaxiangyi@hust.edu.cn，J.Yi@tudelft.nl
clc
clear all
close all
clearvars;
addpath('dace','Infill_criteria','Stop_criteria','Test_Function\test_function_single_fidelity','Failure_Probability_evaluation')
% definition of the teste function four_branches_function
for fun_name={'four_branches_function'}
    test_function=char(fun_name);
    % definition of the learning function
    Infill_critria='Infill_RLCB';
    % import the basic information of the test function
    [num_vari,mu,sigma,design_space,given_thresholds]=test_function_single_fidelity(test_function);
    % number of runs
    for error=given_thresholds
        num_trials=2;
        record.result=zeros(num_trials,5);
        for run=1:num_trials
            fprintf('--------------Test_function:%s; Infill_criteria:%s; run: %d ----------------------------\n',test_function,Infill_critria,run);
            % number of initial sampling
            num_initial_sample=10;
            % design sapce of the constrain optimization problem  define as
            % sampling and calculate the respons
            sample_x = repmat(design_space(1,:),num_initial_sample,1)+repmat(design_space(2,:)-design_space(1,:),num_initial_sample,1).*...
                lhsdesign(num_initial_sample,num_vari,'criterion','maximin','iterations',1000);
            sample_y = feval(test_function,sample_x);
            %  constructing the original Kriging model
            kriging_model = dacefit(sample_x,sample_y,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
            % evaluate the performance of the current iteration
            gen=1;
            iteration=0;
            num_search=10^4*gen;
            search_x=MCS_Population_Generation(mu,sigma,num_search);
            [pf_estimate,pf_real,cov_estimate,real_relative_error]=reliabiliy_evaluation_single_fidelity(search_x,mu,kriging_model,test_function);
            while cov_estimate>0.05 || iteration<=2
                num_search=10^4*gen;
                search_x=MCS_Population_Generation(mu,sigma,num_search);
                stop_value=bootstrap_stop_single_fidelity(search_x,kriging_model,0.05);
                while stop_value>error && iteration<500
                    switch Infill_critria
                        case 'Infill_RLCB'
                            obj=Infill_RLCB(search_x,kriging_model,sample_x);
                            [bestobj,Index]=min(obj);
                            xselected=search_x(Index,:);
                        case 'Infill_H'
                            obj=Infill_H(search_x,kriging_model);
                            [bestobj,Index]=min(obj);
                            xselected=search_x(Index,:);
                        case 'Infill_LIF'
                            obj=Infill_LIF(search_x,kriging_model,mu,sigma);
                            [bestobj,Index]=min(obj);
                            xselected=search_x(Index,:);
                        case 'Infill_eff'
                            z=0;
                            %             [xselected,bestobj]=yichuan(@(x)Infill_eff(x,kriging_model,z),num_vari,design_space(1,:),design_space(2,:));
                            obj=Infill_eff(search_x,kriging_model,z);
                            [bestobj,Index]=min(obj);
                            xselected=search_x(Index,:);
                        case 'Infill_U'
                            %             [xselected,bestobj]=yichuan(@(x)Infill_U(x,kriging_model),num_vari,design_space(1,:),design_space(2,:));
                            obj=Infill_U(search_x,kriging_model);
                            [bestobj,Index]=min(obj);
                            xselected=search_x(Index,:);
                            
                        case'Infill_RU'
                            obj=Infill_RU(search_x,kriging_model,sample_x);
                            [bestobj,Index]=min(obj);
                            xselected=search_x(Index,:);
                    end
                    if ~ismember(sample_x,xselected)
                        fselected=feval(test_function,xselected);
                        sample_x=[sample_x;xselected];
                        sample_y=[sample_y;fselected];
                        kriging_model=dacefit(sample_x,sample_y,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
                    else
                        Mu=repmat(mu,num_search,1);
                        Sigma=repmat(sigma,num_search,1);
                        search_x=normrnd(Mu,Sigma);
                        kriging_model=dacefit(sample_x,sample_y,'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
                    end
                    iteration=iteration+1;
%                     visualization the iteration process
                                if size(sample_x,2)==2
                                    plot_figures_single_fidelity(design_space,test_function,kriging_model,sample_x,num_initial_sample,iteration);
                                end
                    stop_value=bootstrap_stop_single_fidelity(search_x,kriging_model,0.05);
                    % print the current information to the screen
                 [pf_estimate,pf_real,~,real_relative_error]=reliabiliy_evaluation_single_fidelity(search_x,mu,kriging_model,test_function);
                fprintf(' Run=%d; iteration=%d; pf_estimate=%f; pf_real=%f; stop_value=%f  real_error=%f \n', run,iteration ,pf_estimate,pf_real,stop_value,real_relative_error);
                end
                % calclulate the failure probility
                [pf_estimate,pf_real,cov_estimate,real_relative_error]=reliabiliy_evaluation_single_fidelity(search_x,mu,kriging_model,test_function);
                gen=gen+1;
            end
            record.result(run,:)=[pf_estimate,pf_real,cov_estimate,size(sample_x,1),real_relative_error];
        fprintf(' Run=%d; iteration=%d; pf_estimate=%f; pf_real=%f; stop_value=%f  real_error=%f evalution=%d \n', run,iteration ,pf_estimate,pf_real,stop_value,real_relative_error,size(sample_x,1));
            save(strcat('Results/',test_function,'/',mfilename,'_',num2str(100*error),'.mat'),'record');
        end
        record.mean=mean(record.result);
        save(strcat('Results/',test_function,'/',mfilename,'_',num2str(100*error),'.mat'),'record');
        
    end
end
