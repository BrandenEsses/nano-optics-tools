%% Compute the permittivity of SiO2 from Lorentzian fits
%  Input arguments:
%  nu   = optical frequency, in cm^-1
%  T    = temperature, in K
%
%  Output arguments:
%  E    = complex permittivity
%
%  Theoretical model for the T-dependent permittivity:
%
%  E(nu, T) = 1 + F(nu_TO / nu) / nu^2
%
%  where F is some function and nu_TO(T) = f(T) nu_TO(T = 0)
%  is the T-dependent phonon frequency, and
%  f(T) is a linear function descreasing with T.
%  Note that this model preserves the oscillator strength. In reality,
%  it should descrease due to thermal expansion.
%
%  v0.0 11/22/2022 Michael Fogler

function E = SiO2Model4Raschke(nu, T)

%  Phonon softening factor, an adjustable parameter (~3% per 1000 K)
f = 1 - (T - 300) * 3e-5;

nu = nu / f; % rescaled frequency

if true
    load SiO2Model_Re.txt; load SiO2Model_Im.txt;
    
    E = interp1(SiO2Model_Re(:, 1), SiO2Model_Re(:, 2), nu, "spline") ...
      + 1i * interp1(SiO2Model_Im(:, 1), SiO2Model_Im(:, 2), nu, "spline");
else
    E = newSiO2Model(nu);
end

E = 1 + (E - 1) / f^2;
end