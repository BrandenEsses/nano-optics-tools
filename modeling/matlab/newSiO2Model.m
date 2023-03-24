%% Compute the dielectric function of SiO2 from Lorentianz fits
%  precondition:
%    - nu, optical frequency or energy
%  postcondition:
%    - E, complex dynamic dielectric function
function Enew = newSiO2Model(nu)
nu = real(nu); % ignore imaginary part.
% % % lp = {[sqrt(14.258 * 1072 * 47.727), 1072, 47.727], ...
% % %     [sqrt(6.0733 * 1046 * 4.1927), 1046, 4.1927], ...
% % %     [sqrt(0.19118 * 1181 * 22.151), 1181, 22.151], ...
% % %     [sqrt(-0.16911 * 1270 * 206.35), 1270, 206.35], ...
% % %     [sqrt(0.42108 * 1210 * 74.021), 1210, 74.021]};
% % % eps0 = 1.9259;
% % % 
% % % E = ones(size(nu)) * eps0;
% % % for z = 1:numel(lp)
% % %     E = E + lp{z}(1)^2 ./ (lp{z}(2)^2 - nu.^2 - 1i * lp{z}(3) * nu);
% % % end

% is working but without 800cm-1 phonon
% lp = [ 14.261   1072    49.444
%      -0.18416  1270    215.53
%       0.46239  1205    77.698];
% lp(:,3) = lp(:,3) * 0.9;
% eps0 = 1.9259;
% E = ones(size(nu)) * eps0;

% for z = 1: size(lp, 1)
%     E = E + lp(z, 1) * lp(z, 2) * lp(z, 3) ...
%       ./ (lp(z, 2)^2 - nu.^2 - 1i * lp(z, 3) * nu);
% end



% new dispersion inlcuding 800cm-1 phonon from Kucirkova+Navratil (1994) 
% with Voigt-like profile according to paper; in part translated for 
% matlab from Alex model

lp2=[1.73  1147  40  81   %amplitude position Lorentz- Gauss-width [cm-1]
     1.5   1091  .3  24
     3.3   1060   8  30
      .35   808  10  32]; 

%other values from same paper as alternative:
%lp2=[1.7 1150 40 80
%     3.1 1080 .1 26
%     1.1 1045 10 17
%     .34  800 .1 56];  
  
%compute correct prefactors:
lp2(:,1) = lp2(:,1)*1e5*sqrt(pi)/2;
lp2(:,1) = lp2(:,1)./lp2(:,2); %divide amplitude by frequency
lp2(:,1) = lp2(:,1)./lp2(:,4); %devide amplitude by Gauss-width

eps0new = 1.96;
Enew = ones(size(nu)) * eps0new;

for n=1: size(lp2,1)
    wminus=faddeeva((nu-lp2(n,2))/(lp2(n,4))+1i*lp2(n,3)/(2*lp2(n,4)));
    wplus=faddeeva((nu+lp2(n,2))/(lp2(n,4))+1i*lp2(n,3)/(2*lp2(n,4)));
    Enew=Enew + lp2(n,1)*1i*(wminus-wplus);
end

% plot(nu,real(E),nu,real(Enew),nu,imag(E),nu,imag(Enew));
% savefile = ['imagPartOfOLDMishaModel'];
% savematrix = [nu' imag(E)'];
% save(savefile, '-ASCII', 'savematrix');
end


% computes the complex probability function representing the Voigt profiles
% according to Kucirkova & Navratil (1994) and equivalent to the Faddeeva
% function w(z)=exp(-z^2)*erfc(-iz) for imag(z)>0
% see Weideman "computation of the complex error function" (1994) for
% Matlab code
function result = faddeeva(z)
    N=50; %rational series of N terms
    M=2*N; 
    M2=2*M;
    k=[-M+1:1:M-1]';
    L=sqrt(N/sqrt(2));
    theta=k*pi/M;
    t=L*tan(theta/2);
    f=exp(-t.^2).*(L^2+t.^2);
    f=[0; f];
    a=real(fft(fftshift(f)))/M2;
    a=flipud(a(2:N+1));
    Z=(L+1i*z)./(L-1i*z);
    p=polyval(a,Z);
    w=2*p./(L-1i*z).^2+(1/sqrt(pi))./(L-1i*z);
    result=w;
end
