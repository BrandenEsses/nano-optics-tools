%% function sNbar = Spheroid4Raschke(nu, T, Nd)
%
%  Superstrate (air)
%  Drude metal, thickness d_1
%  SiO2, thickness d_2
%  Substrate (Si)
%
%  Input arguments:
%
%  nu       = frequency in cm^-1, a vector of length Nnu
%  Nd       = demodulation order, e.g., 3 or 4
%  T        = temperature, in K, a row vector of length NT
%             The 1st elelement should be the reference T = 300 K
%
%  Output arguments:
%
%  sNbar    = normalized and demodulated complex SNOM signal, 1 x Nnu
%
%  v0.0 10/03/2022 Michael Fogler
%  v0.1 11/22/2022 changed to varying T instead of omega_p, gamma_D
%
%  Example of usage:
%  T = [300 800]; nu = linspace(800,1600,401); Spheroid4Raschke(nu,T,2);

function [sNbar,amplitude,phase] = Spheroid4Raschke(nu, T, Nd)
%% Parameters
cm = 1; nm = 1e-7 * cm; % length units

%  Geometry parameters:
%  These parameters are adjustable. However, don't go crazy!
%  Making big changes may give meaningless results. 
%
a = 35 * nm;    % radius of curvature of tip apex
L = 10 * a;     % long semi-axis of the spheroid, L > a
dz = 40 * nm;   % tip oscillation amplitude


%  Sample parameters:
%  These parameters are adjustable unless they are precisely known
%  from experiment
d1 = 40 * nm;   % thickness of medium 1 (Drude metal)
d2 = 300 * nm;  % thickness of medium 2 (SiO2)
omega_p  = 1001;% bare plasma frequency, in cm^-1, a vector of length Nmu
gamma_D  = 1000;% broadening parameter of the Drude conductivity, in cm^-1

%%
%% The following part is the core algorith that needs no changes
%% The "Make quick plots" section can be modified as desired
%%

%  Spheroid geometry parameters
F = sqrt(L * (L - a));  % focal length of the spheroid
xi = L / F;             % auxiliary variable

%  Tip-sample separation, a row vector
%  Use Gauss-Chebyshev quadrature nodes and weights (a column vector)
%  to implement the demodulation procedure
Nz = 41;
z_tip_min = 0.02 * a;                   % min tip-sample distance
x = pi / 2 / Nz * (1: 2: 2 * Nz - 1);   % equidistant "angle" grid
z = z_tip_min + dz * (1 - cos(x));      % z-positions of the tip
% Wdem = cos(double(Nd) * x).' / double(Nz);              % (weights) * (demod factor)
Wdem = cos(Nd * x).' / Nz;
%% Compute matrix elements
% tic;
% fprintf('\n Calculating matrix elements\n')
Nm = 10; % dimension of the matrix H

%  Matrix element H_{ln} = \int f(q) dq from 0 to infinity, where
%  x = q * F,
%
%  f(q) = 2 * pi * I_(l + 1/2)(q F) * I_(n + 1/2)(q F) / q
%       * r_P(q) * exp(-2 q (z + L - F))),
%
%  and I_a(z) = besseli(a, z, 1) = Mathematica's BesselI[a, z] Exp[-z].
%
%  Here we calculate only the common factor
%
%  H0_{ln} = 2 * pi * I_(l + 1/2)(q F) * I_(n + 1/2)(q F) / q
%
%  The remaining two factors (the exponential and the r_P)
%  are added later in function Chi0Spheroid()
%
%  The integral is approximated by the sum(f(q_j) * W_j)
%  where q_j are the Gauss-Legendre nodes, W_j are the weights

%  Compute the momentum grid for the Gauss-Legendre quadrature
Nq = 500;       % no. nodes of the quadrature
q_max = 16 / a; % upper limit on q in the integral (not inf, of course)

%  Construct the Legendre polynomials matrix
tmp = (1: Nq - 1);
diagnl = tmp ./ sqrt(4 * tmp.^2 - 1);
Bin = diag(diagnl, 1) + diag(diagnl, -1);

%  Diagonalize it and sort the eigensystem
[P, M] = eig(Bin);
tmp = diag(M); % roots of P_{Nq}
[q, idx] = sort(tmp); W = P(1, idx).';

%  Map from (-1, 1) to (0, q_max)
q = q_max / 2 * (q + 1);    % Gauss-Legendre nodes
W = q_max * W.^2;           % Gauss-Legendre weights

%  Proceed to calculating matrix elements H01, H02, H03
qF = q * F;
H01 = zeros(Nq, Nm);
for n = 1: Nm
        nn = n + 0.5;
        ll = n + 0.5;
        
        B = besseli(nn, qF, 1) .* besseli(ll, qF, 1);
        H01(:, n) = 2 * pi ./ q .* B .* W;    
end
H02 = zeros(Nq, Nm - 1);
for n = 1: Nm - 1
        nn = n + 0.5;
        ll = n + 1.5;
        
        B = besseli(nn, qF, 1) .* besseli(ll, qF, 1);
        H02(:, n) = 2 * pi ./ q .* B .* W;    
end
H03 = zeros(Nq, Nm - 2);
for n = 1: Nm - 2
        nn = n + 0.5;
        ll = Nm + 0.5;
        
        B = besseli(nn, qF, 1) .* besseli(ll, qF, 1);
        H03(:, n) = 2 * pi ./ q .* B .* W;    
end

%  Compute Lambda (vector of the spheroid's polarizabilities)
%
%  First, compute LegendreQ ratios QQ_n = Q_{n + 1} / Q_n
%  Starting value of QQ_n is the asymptotic limit, which we impose at
%  high enough n = Nm + n_max. Lower qq_n's are by backward recursion
%
%  Transient values: from n = Nm + n_max to n = Nm - 1
tmp = 1 / (xi + sqrt(xi^2 - 1));
n_max = -Nm / (0.01 + log(tmp));
for n = Nm + n_max - 1: -1: Nm - 1
    tmp = (n + 1) / ((2 * n + 3) * xi - (n + 2) * tmp);
end

%  Useful values: from n = 1 to Nm - 2
QQ = zeros(1, Nm - 1); QQ(Nm - 1) = tmp;
for n = Nm - 2: -1: 1
    QQ(n) = (n + 1) / ((2 * n + 3) * xi - (n + 2) * QQ(n + 1));
end

%  Compute LegendreP ratios PP_n = P_{n} / P_{n + 1}.
%  Starting value is PP_0 = 1 / xi. Higher pp_n's are by forward recursion
PP = zeros(1, Nm - 1); PP(1) = 2 / (3 * xi - 1 / xi);
for n = 2: Nm - 1
    PP(n) = (n + 1) / ((2 * n + 1) * xi - n * PP(n - 1));
end

%  The cumulative product of (QQ * PP) gives Lambda (up to 4 /(2 n + 1) )
Q1_P1 = 0.5 * log((xi + 1) / (xi - 1)) - 1 / xi;
Lambda = cumprod([Q1_P1, QQ .* PP]);
Lambda = Lambda * 4 ./ (3: 2: 2 * Nm + 1);
% A side note: dipole polarizability is
% p_inf = F^3 / (2.25 * Lambda(1));

irtLambda = 1 ./ sqrt(Lambda); % is used for rescaling matrix H below
% toc;

%% Compute the tip polarizability chi0
fprintf('\n Calculating polarizability chi0\n')
tic;

NT = length(T);
Nnu = length(nu);
chi0 = zeros(Nnu, Nz, NT);

for j = 1: NT       
    % rP4Raschke() returns Nnu x Nq array
    rp = rP4Raschke(q, nu, omega_p, gamma_D, T(j),...
             d1, d2);

    % Chi0Spheroid() takes Nnu x Nq array and returns
    % a Nnu x Nz array whose 1st row corresponds to z = z_tip_min
    chi0(:, :, j) = Chi0Spheroid(z, a, L, nu, rp,...
                      H01, H02, H03, irtLambda, q);
end

%  Reference material (Au)
rpref = 0.97 * ones(1, Nq); % r_P of the reference material at nu(1)

chi0ref = Chi0Spheroid(z, a, L, nu(1), rpref,...
          H01, H02, H03, irtLambda, q); % chi0 of the reference material
toc;

%% Do the demodulation
% tic;
% fprintf('\n Demodulating\n')

%  Demodulate the SNOM signal for the sample
%  Use the Gauss-Chebyshev quadrature
chiN = zeros(Nnu, NT);
for j = 1: NT
        chiN(:, j) = chi0(:, :, j) * Wdem;
end
%  Do the same for the reference material
chiNref = chi0ref * Wdem;
% toc;

%% Add far-field factors, assuming angle of incidence theta = pi / 4
%  Note: FF momentum is (omega/c) sin(theta) = 2 pi nu / sqrt(2)
sN = zeros(Nnu, NT);
for j = 1: NT
    rpFF = rP4Raschke(2 * pi * nu / sqrt(2), nu, omega_p, gamma_D, ...
           T(j), d1, d2);   
    rpFF = diag(rpFF); % a lazy way to get requisite (q, nu) pairs
    FFF = (1 + rpFF).^2;

    sN(:, j) = chiN(:, j) .* FFF;
end

%  FFF and SNOM signal for the reference material
rpFFref = rpref(1, 1);
FFFref = (1 + rpFFref).^2;
sNref = chiNref * FFFref;

sNbar = sN / sNref;

%% Prepare phase for plotting
phase_0 = phase_unwrap(angle(sNbar(:, 1)));
phase = phase_unwrap(angle(sNbar));
phase_dif = phase_unwrap(angle(sNbar(:, 2)) - angle(sNbar(:, 1)));

amplitude = abs(sNbar);

end

%% Phase unwrapping
%  This auxiliary function makes the phase continuous
%  (eliminates 2\pi-jumps) to make plots look better
function p = phase_unwrap(phi)
N = numel(phi);
p = phi;
for j = 1: N - 1
    dp = p(j + 1) - p(j);
    if abs(dp + 2 * pi) < abs(dp)
        p(j + 1: N) = p(j + 1: N) + 2 * pi;
    elseif abs(dp - 2 * pi) < abs(dp)
        p(j + 1: N) = p(j + 1: N) - 2 * pi;
    end
end
p = p / pi * 180;
end