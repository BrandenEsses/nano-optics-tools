%% function rP = rP4Raschke(q, nu, omega_p, gamma_D, T, d1, d2)
%
%  Calculates the P-polarization refectivity r_P of the following
%  multi-layer system:
%
%  Superstrate (air)
%  Drude metal, thickness d_1
%  SiO2, thickness d_2
%  Substrate (Si)
%
%  Input arguments: 
%
%  q        = momenum, length Nq vector
%  nu       = frequency, in cm^-1, length Nnu vector
%  omega_p  = plasma frequency, in cm^-1 
%  gamma_D  = broadening parameter of the Drude conductivity
%  T        = temperature, in K
%  d1       = thickness of the Drude metal
%  d2       = thickness of SiO2
%
%  Output arguments:
%
%  rP       = P-polarization refectivity, Nnu x Nq array
%
%  v0.0 10/03/2022 Michael Fogler
%  v0.1 11/22/2022 Michael Fogler 

function rP = rP4Raschke(q, nu, omega_p, gamma_D, T, d1, d2)

%  Parameters
eps_inf = 13;       % high-frequency permittivity of the Drude metal 
gamma_inf = 1e-6;   % infinitesimal damping

%  Grid of (q, omega)
%  Note the 2 * pi, the units conversion from (frequency) / (2 pi c),
%  measured in cm^-1, to (frequency) / c.
[q_v, w_v] = meshgrid(q, 2 * pi * nu + gamma_inf * 1i);

%% Individual layers and interfaces
% Layer 0: air
k0 = sqrt(w_v.^2 - q_v.^2);
% idx = imag(k0) < 0; k0(idx) = -k0(idx);
Q0 = 1 ./ k0;

% Layers 1: Drude metal
E1 = eps_inf - omega_p^2 ./ (nu + 1i * gamma_D) ./ nu;
[~, E1] = meshgrid(q, E1);
k1 = sqrt(E1 .* w_v.^2 - q_v.^2); % k^z
% idx = imag(k1) < 0; k1(idx) = -k1(idx);
Q1 = E1 ./ k1;

% Layer 2: SiO2
E2 = SiO2Model4Raschke(nu, T); % dielectric function of SiO2
[~, E2] = meshgrid(q, E2);
k2 = sqrt(E2 .* w_v.^2 - q_v.^2);
% idx = imag(k2) < 0; k2(idx) = -k2(idx);
Q2 = E2 ./ k2;

% Layer 3: substrate (Si)
E3 = 11.7; % dielectric function of the substrate (Si)
[~, E3] = meshgrid(q, E3);
k3 = sqrt(E3 .* w_v.^2 - q_v.^2);
% idx = imag(k3) < 0; k3(idx) = -k3(idx);
Q3 = E3 ./ k3;

% Compute 23 reflection coefficient
r23 = (Q3 - Q2) ./ (Q3 + Q2);

% Compute 12 reflection coefficients
r12 = (Q2 - Q1) ./ (Q2 + Q1);
r21 = -r12;

% Compute 01 reflection coefficients
r01 = (Q1 - Q0) ./ (Q1 + Q0);
r10 = -r01;

%% Implement the multi-layer recurrence relation
% Starting value: the bottom (23) interface
rP = r23;

% Step:
% Add the phase factor due to roundtrip distance
rP = rP .* exp(2i * d2 * k1); % note: k2 = k1
% Do the fraction
% R = (r12 + (1 - r12 - r21) .* R) ./ (1 - r21 .* R);
rP = (r12 + rP) ./ (1 - r21 .* rP); % simpler form

% Step:
% Add the phase factor due to roundtrip distance
rP = rP .* exp(2i * d1 * k1);
% Do the fraction
% R = (r01 + (1 - r01 - r10) .* R) ./ (1 - r10 .* R);
rP = (r01 + rP) ./ (1 - r10 .* rP); % simpler form

% We are done
end