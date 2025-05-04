# BQ_IRKA
"This repository contains MATLAB code for numerical experiments performed in the thesis "H2-optimal model reduction of bilinear systems with quadratic output" by Octavius Guenther at TU Berlin.

## requirements
-only MATLAB R2024b

## installation
- copy the repository
- TO OBTAIN ALL THE PLOTS IN THE THESIS run ./code/run_all.m in MATLAB 2024b


## repository structure
- './code/' MATLAB source code
- './plots/' generated plots

## tutorial for playing around with BQ-systems
- create a BQ-system via the console via with "Sigma = BQ_system(dim,model_name)", where modelname is either 'heat','flow','rand'.
- run "Sigma_opt=BQ_IRKA_v3(Sigma,r)", where r is the size of desired reduced-order system
- to solve the generalized Sylvester or Lyapunov-equations iteratively run "gen_sylv_naive(Sigma,Sigma,max_iter,tol)"

