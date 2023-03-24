%% Chi0Spheroid.m
%  function chi0 = Chi0Spheroid(z, a, L, nu, r_P, H01, H02, H03,...
%                  irtLambda, q)
%  Computes the total dipole moment for a grid of (nu, z) points
%
%  Input arguments:
%
%  z        = tip-sample separation/radius, a vector of length Nz
%  a        = tip radius
%  L        = long spheroid semi-axis/radius
%  w        = frequency (cm^-1), a vector of length Nnu
%  r_P      = matrix of precalculated r_P, size Nnu x Nq
%  H01      = diagonal matrix elements, Nq x Nm
%  H02      = matrix elements of 1st upper diagonal, Nq x (Nm - 1)
%  H03      = matrix elements of the last column, Nq x (Nm - 2)
%  irtLambda= 1 ./ sqrt(Lambda) is used for rescaling the matrix H
%  q        = momenta, a column vector of length Nq
%
%  Output argument:
%
%  chi0     = dipole moment, Nnu x Nz array
%
%  v0.0 10/02/2022 Michael Fogler

function chi0 = Chi0Spheroid(z, a, L, nu, r_P, H01, H02, H03,...
                irtLambda, q)
%% Spheroid geometry parameters
F = sqrt(L^2 - L * a);  % focal distance

%% 
Nz = length(z);
Nnu = length(nu);
Nm  = size(H01, 2); % dimension of matrix H
chi0 = zeros(Nnu, Nz);

z_minus_F = z + L - F; % the F-subtraction is because of besseli(..., 1)
        
for i = 1: Nnu
    beta = r_P(i, :).'; % column vector of length Nq
    H01_ = H01 .* beta; % old: (beta * ones(1,Nm));
    H02_ = H02 .* beta; % old: (beta * ones(1,Nm-1));
    H03_ = H03 .* beta; % old: (beta * ones(1,Nm-2));
    
    for j = 1: Nz
        exp_fac = exp(-2 * z_minus_F(j) * q); % column vector of length Nq
        H01__ = H01_ .* exp_fac; % old: (exp_fac * ones(1,Nm));
        H02__ = H02_ .* exp_fac; % old: (exp_fac * ones(1,Nm-1));
        H03__ = H03_ .* exp_fac; % old: (exp_fac * ones(1,Nm-2));

        % Do the quadrature over q by summing over the rows
        HH1 = sum(H01__); 
        HH2 = sum(H02__);
        HH3 = sum(H03__);
        
        %  main diagonal
        H = diag(HH1);

        %  1st diagonal
        for n = 1: Nm - 1
            H(n, n + 1) = HH2(n);
        end

        %  last column
        for n = 1: Nm - 2
            H(n, Nm) = HH3(n);
        end

        %  the rest of the upper triangle is found by recursion
        for g = Nm - 2: -1: 2                  
            for k = 1: g - 1
                n  = g - k + 1;
                n2 = Nm - k;
                H(n - 1, n2) = H(n + 1, n2)...
                             - (2 * n + 1) / (2 * n2 + 1)...
                             * (H(n, n2 + 1) - H(n, n2 - 1));
            end
        end

        %  Scale H by 1 / sqrt(Lambda) to lower the Rcond number
        for n = 1: Nm
            for k = n: Nm        
                H(n, k) = H(n, k) * irtLambda(n) * irtLambda(k);
                H(k, n) = H(n, k);
            end
        end

        %  Solve for the dipole moment
        M = eye(Nm) - H;

%         warning('off', 'MATLAB:nearlySingularMatrix');
        A = M \ [4 / 3 * irtLambda(1)^2; zeros(Nm - 1, 1)];
%         warning('on', 'MATLAB:nearlySingularMatrix');

        chi0(i, j) = A(1) * ((F/a)^3 / 3);
    end
end
end