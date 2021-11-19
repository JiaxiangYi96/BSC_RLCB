function [num_vari,mu,sigma,design_space,stopping_thresholds] = test_function_single_fidelity(name)
%--------------------------------------------------------------------------
% single-objective optimization problem
switch name
    
    case 'four_branches_function'
        num_vari=2;   mu=[0,0];sigma=[1,1];design_space=[mu-5*sigma;mu+5*sigma]; stopping_thresholds=[0.03 0.02 0.01];
    case 'four_branches_function_2'
        num_vari=2;   mu=[0,0];sigma=[1,1];design_space=[mu-5*sigma;mu+5*sigma];stopping_thresholds=[0.03 0.02 0.01];
     case 'multimodal_function'
        num_vari=2;   mu=[1.5,2.5];sigma=[1,1];design_space=[mu-5*sigma;mu+5*sigma];stopping_thresholds=[0.03 0.02 0.01];
      case 'nonlinear_oscillator_function'
        num_vari=6;   mu=[1,1,0.1,0.5,1,1];sigma=[0.05,0.1,0.01,0.05,0.1,0.1];design_space=[mu-5*sigma;mu+5*sigma];  stopping_thresholds=[0.03 0.02 0.01];
       case 'roof_truss'
          num_vari=6;  
          mu=[20000,13.5,9.82*10^-4,0.04,1.2*10^11,3*10^10];
          sigma=[2000,0.50,6*10^-5,0.008,1.2*10^10,3*10^9];
          design_space=[mu-5*sigma;mu+5*sigma]; 
          stopping_thresholds=[0.03 0.02 0.01];
    case'highly_nonlinear_log'
        num_vari=6;   
        mu=[120 120 120 120 50 40];
        sigma=[12 12 12 12 15 12];  
       design_space=[mu-5*sigma;mu+5*sigma];  
       stopping_thresholds=[0.03 0.02 0.01];
    otherwise
        sprintf('%s function is not defined', name);
end

end
