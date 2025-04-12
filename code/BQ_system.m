classdef BQ_system 
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        A
        N
        b
        M
        dim
    end
                                                 
    methods
        function obj = BQ_system(n,model,varargin)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            switch model

                case 0

                    obj.A = sparse(n,n);
                    obj.N = sparse(n,n);
                    obj.b = sparse(n,1);
                    obj.M = sparse(n,n);
                    obj.dim = n;

                case 'heat'
                    [obj.A,obj.N,obj.b,obj.M] = heat_equation_discretization_single(n);
                    obj.dim = size(obj.A,1);
                    

                case 'rand'
                    obj.A = randi([-10 10],n,n);
                    obj.N = randi([-10 10],n,n);
                    obj.b = randi([-10 10],n,1);
                    M = randi([-10 10],n,n);
                    
                    obj.M = M'*M;
                    obj.dim = size(obj.A,1);
                    
                    %make it stable, i.e. eig(A) < 0 
                    obj.A = obj.A - 2*abs(eigs(obj.A,1))*eye(obj.dim);

                case 'flow'
                    [obj.A,obj.N,obj.b,obj.M] = flow_model_discretization(n);
                    obj.dim = size(obj.A,1);
            end
            
        end
        
        function b = stability(obj)
            p = poles(obj,'all'); 
            b = all(p(:)<0);
        end

        % returns the the poles of the system
        function p = poles(obj,all)
            if (~exist('all','var'))
                p = eig(full(obj.A));
            else
            p = eigs(obj.A);
            end
        end
        
        %in case of convergence the optimal interpolation points are given
        %as -sigma(A_r)
        function p = optimal_points(obj)
            p = eig(full(-1*obj.A));
        end
        
        %calculate output y
        %creates a pair [t,Y] that can be plotted
        %initial value x0=0
        %solver used ode45
        function [t,Y] = get_output(obj,fun_input,grid)
            %dimensio of the object
            n = obj.dim;
            

            %function handle for the solver
            f = @(t,x) obj.A*x + obj.N*x*fun_input(t) + obj.b*fun_input(t);

            %initial value
            x0 = zeros(n,1);

            %solve in [0,2]
            [t,X] = ode45(f,grid,x0);

            %size of the discretization
            N = size(t,1);
            %calculating output Y
            Y = arrayfun(@(i) X(i,:)*obj.M*X(i,:)',1:N);

            %Y has to be transposed 
            [t,Y'];
        end

        %Gramian convergence condition
        function  a = gramian_conv(obj)
            a = sqrt(norm(full(obj.N*obj.N')));
        end

        %computes the Gramians by solving the corresponsing Lyapunov
        %euqations
        function [P,Q] = getGramians(obj,method)
            if (~exist('method','var'))
                method = 'direct';
            end
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            P = gen_sylv(obj,obj,method);
            Q = gen_sylvZ(obj,obj,P,method);
        end
        
        %computes only the reachability Gramian
        function P = getReachGram(obj,method)
            if (~exist('method','var') )
                method = 'direct';
            end
                P = gen_sylv(obj,obj,method);
        end
        
        %computes the norm |Sigma|_H2 of a BQ system
        %opt denotes one of the two formulac that is used 
        %gram denotes the full order Gramian P that can be passed to omit
        %calculating large Lyapunov equations
        function norm = getH2norm(obj,opt,gram,method)
            if (~exist('method','var'))
                method = 'direct';
            end

            if (~exist('opt', 'var'))
               opt = 'P';
            end

            switch opt
                case 'P'

                %norm = trace(MPMP)
                if (~exist('gram','var'))

                    P = getReachGram(obj,method);
                else
                    P = gram;
                end
                norm = (P*obj.M)^2;
                norm = trace(norm);

                case 'Q'
                
                %norm = b'*Q*b      
                [~,Q] = getGramians(obj,method);
                norm = full(obj.b'*Q*obj.b);

            end


        end


        %computes the error | Sigma1 - Sigma2 |
        %opt denotes the method which formula to apply for calculating the
        %norm
        % FOM_gram can be passed for not computing the large full order
        % gramian iteratively
        function error = getErrorH2norm(obj1,obj2,opt,FOM_gram)

            if (~exist('opt', 'var'))
                    opt = 'P';
            end

            switch opt
                case 'P'
                
                if(~exist('FOM_gram','var'))

                    P1 = getReachGram(obj1);
                else
                    P1 = FOM_gram;
                end
                
                P2 = getReachGram(obj2);
                X = gen_sylv(obj1,obj2);
    
                error =trace(P1*obj1.M*P1*obj1.M)-2*trace(X*obj2.M*X'*obj1.M) + trace(P2*obj2.M*P2*obj2.M);

                case 'Q'

                [~,Q1]=getGramians(obj1);
                [~,Q2]=getGramians(obj2);
            
                X = gen_sylv(obj1,obj2);
                Z = gen_sylvZ(obj1,obj2,-1*X);

                error = full(obj1.b'*Q1*obj1.b + 2*obj1.b'*Z*obj2.b + obj2.b'*Q2*obj2.b);

            end
            
        
        end

        % computes the spectral radius of the 
        % generalized sylvester operator LA^-1*Pi
        % corresponding to the sylvester equation 
        % AX + XA_r^T + NXN_r^T = - bb_r^T
        % If the spectral radius is < 1, the naive iteration for solving
        % generalized Sylvester equations will converge
        
        function a = spec_rad(Sigma_1,Sigma_2)
        n = Sigma_1.dim;
        r = Sigma_2.dim;
        
        LA = kron(speye(r),Sigma_1.A)+kron(Sigma_2.A,speye(n));
        
        Pi = sparse(kron(Sigma_2.N,Sigma_1.N));
        
        
        C = LA\Pi;
        
        a = abs(eigs(C,1));

        end


        %displays the BQ system py printing all system matrices
        function disp(obj)
            fprintf('A')
            disp(full(obj.A))
            fprintf('N')
            disp(full(obj.N))
            fprintf('b')
            disp(full(obj.b))
            fprintf('M')
            disp(full(obj.M))
        end

        function size(obj)
            size(obj.A)
        end
        
    end
end