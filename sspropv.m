function [u1x,u1y] = sspropv(u0x,u0y,dt,dz,nz,alphaa,alphab,betapa,betapb,gamma,psp,method,maxiter,tol)
% 该函数使用分步傅立叶方法求解光纤中脉冲传播的耦合模式非线性薛定谔方程。
% 模型中包括以下影响：群速度色散(GVD)/高阶色散/偏振相关损耗(PDL)/任意光纤双折射/自相位调制(SPM)
% USAGE
%
% [u1x,u1y] = sspropv(u0x,u0y,dt,dz,nz,alphaa,alphab,betapa,betapb,gamma);
% [u1x,u1y] = sspropv(u0x,u0y,dt,dz,nz,alphaa,alphab,betapa,betapb,gamma,psp);
% [u1x,u1y] = sspropv(u0x,u0y,dt,dz,nz,alphaa,alphab,betapa,betapb,gamma,psp,method);
% [u1x,u1y] = sspropv(u0x,u0y,dt,dz,nz,alphaa,alphab,betapa,betapb,gamma,psp,method,maxiter;
% [u1x,u1y] = sspropv(u0x,u0y,dt,dz,nz,alphaa,alphab,betapa,betapb,gamma,psp,method,maxiter,tol);
%
% INPUT
% u0x, u0y        Starting field amplitude components
% dt              Time step
% dz              Propagation step size
% nz              Number of steps to take (i.e. L = dz*nz)
% alphaa, alphab  Power loss coefficients for the two eigenstates(see note (2) below)         
% betapa, betapb  Dispersion polynomial coefs, [beta_0 ... beta_m] for the two eigenstates (see note (3) below)
% gamma           Nonlinearity coefficient
% psp             Polarization eigenstate (PSP) of fiber, see (4)
% method          Which method to use, either 抍ircular?or 抏lliptical? (default = 抏lliptical? see instructions)
% maxiter         Max number of iterations per step (default = 4)
% tol             Convergence tolerance (default = 1e-5)
% 
% OUTPUT
% u1x, u1y        Output field amplitudes
% 
% psp：光纤的偏振本征态。如果psp是标量，给出线性双折射轴相对于x-y轴的方向。
% 如果psp是长度为2的矢量，即psp=[psi，chi]，则它描述椭圆取向和第一偏振本征态的椭圆率。
% 具体来说，(2*psi，2*chi)是庞加莱球面上主本征态的经度和纬度


if (nargin<10)
    error('Not enough input arguments')
end

if (nargin<11)
    psp = [0,0];
end

if (nargin<12)
    method='circular';
end

if (nargin<13)
    maxiter = 4;
end

if (nargin<14)
    tol = 1e-5;
end

nt = length(u0x);
w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]/(dt*nt) ;

if isscalar(psp)
    psi = psp(1);               % Orientation of birefringent axes
    chi = 0;                    % (linear birefringence)
else
    psi = psp(1);               % Orientation of polarization ellipse
    chi = psp(2);               % Ellipticity parameter
end
%% 
if (length(alphaa) == nt)     % If the user manually specifies alpha(w)
    ha = -alphaa/2;
else
    ha = 0;
    for ii = 0:length(alphaa)-1
        ha = ha - alphaa(ii+1)*(w).^ii/factorial(ii);
    end
    ha = ha/2;
end

if (length(betapa) == nt)     % If the user manually specifies beta(w)
    ha = ha - 1j*betapa;
else
    for ii = 0:length(betapa)-1
        ha = ha - 1j*betapa(ii+1)*(w).^ii/factorial(ii);
    end
end

ha = exp(ha.*dz/2); % ha = exp[(-alphaa/2 - j*betaa)*dz/2])

if (length(alphab) == nt)
    hb = -alphab/2;
else
    hb = 0;
    for ii = 0:length(alphab)-1
        hb = hb - alphab(ii+1)*(w).^ii/factorial(ii);
    end
    hb = hb/2;
end

if (length(betapb) == nt)
    hb = hb - 1j*betapb;
else
    for ii = 0:length(betapb)-1
        hb = hb - 1j*betapb(ii+1)*(w).^ii/factorial(ii);
    end
end

hb = exp(hb.*dz/2);  % hb = exp[(-alphab/2 - j*betab)*dz/2])
%% 圆
if strcmp(method,'circular')   %% CIRCULAR BASIS METHOD %%
    
    % First, rotate coordinates to circular basis:
    u0a = (1/sqrt(2)).*(u0x + 1j*u0y);
    u0b = (1/sqrt(2)).*(1j*u0x + u0y);
    
    % Propagation matrix for linear calcuations
    
    h11 = ( (1+sin(2*chi))*ha + (1-sin(2*chi))*hb )/2;
    h12 = -1j*exp(+1j*2*psi)*cos(2*chi)*(ha-hb)/2;
    h21 = +1j*exp(-1j*2*psi)*cos(2*chi)*(ha-hb)/2;
    h22 = ( (1-sin(2*chi))*ha + (1+sin(2*chi))*hb )/2;
    
    u1a = u0a;
    u1b = u0b;
    uafft = fft(u0a);
    ubfft = fft(u0b);

    for iz = 1:nz
    % Calculate 1st linear half
        uahalf = ifft(h11.*uafft + h12.*ubfft);
        ubhalf = ifft(h21.*uafft + h22.*ubfft);
        for ii = 1:maxiter
          % Calculate nonlinear section
            uva = uahalf .* exp((-1j*(1/3)*gamma*dz).* ...
                                (2*(abs(u0a).^2+abs(u1a).^2)/2 + ...
                                 4*(abs(u0b).^2+abs(u1b).^2)/2));
            uvb = ubhalf .* exp((-1j*(1/3)*gamma*dz).* ...
                                 (2*(abs(u0b).^2+abs(u1b).^2)/2 + ...
                                  4*(abs(u0a).^2+abs(u1a).^2)/2));
            uva = fft(uva);
            uvb = fft(uvb); 
            % Calculate 2nd linear half
            uafft = h11.*uva + h12.*uvb;
            ubfft = h21.*uva + h22.*uvb;
            uva = ifft(uafft);
            uvb = ifft(ubfft);
          
            if ((sqrt(norm(uva-u1a,2).^2+norm(uvb-u1b,2).^2) / ...
               sqrt(norm(u1a,2).^2+norm(u1b,2).^2))<tol)
            % tolerances met, break loop
                u1a = uva;
                u1b = uvb;
                break;
            else
            % tolerances not met, repeat loop
                u1a = uva;
                u1b = uvb;
            end
        end %end convergence iteration
        if (ii == maxiter)
            warning(sprintf('Failed to converge to %f in %d iterations', tol,maxiter));
        end
        u0a = u1a;
        u0b = u1b;
    end %end step iteration
  
  % Rotate back to x-y basis:
    u1x = (1/sqrt(2)).*(u1a-1j*u1b) ;
    u1y = (1/sqrt(2)).*(-1j*u1a+u1b) ;
  
%% 椭圆
elseif strcmp(method,'elliptical')    %% ELLIPTICAL BASIS METHOD %%
  % First, rotate coordinates to elliptical basis of eigenstates:

    u0a = (cos(psi)*cos(chi) - 1j*sin(psi)*sin(chi))*u0x + ...
          (sin(psi)*cos(chi) + 1j*cos(psi)*sin(chi))*u0y;
    u0b =(-sin(psi)*cos(chi) + 1j*cos(psi)*sin(chi))*u0x + ...
          (cos(psi)*cos(chi) + 1j*sin(psi)*sin(chi))*u0y;
  
    u1a = u0a;
    u1b = u0b;
    uafft = fft(u0a);
    ubfft = fft(u0b);

    for iz = 1:nz
    % Calculate 1st linear half
        uahalf = ifft( ha.*uafft );
        ubhalf = ifft( hb.*ubfft );
        for ii = 1:maxiter
          % Calculate nonlinear section
            uva = uahalf.*exp((-1j*(1/3)*gamma*dz).* ...
                              ((2+cos(2*chi)^2)*(abs(u0a).^2+abs(u1a).^2)/2 + ... 3
                               (2+2*sin(2*chi)^2)*(abs(u0b).^2+abs(u1b).^2)/2));
            uvb = ubhalf.*exp((-1j*(1/3)*gamma*dz).* ...
                              ((2+cos(2*chi)^2)*(abs(u0b).^2+abs(u1b).^2)/2 + ...
                               (2+2*sin(2*chi)^2)*(abs(u0a).^2+abs(u1a).^2)/2));
            uva = fft(uva);
            uvb = fft(uvb);
            % Calculate 2nd linear half
            uafft = ha.*uva;
            ubfft = hb.*uvb;
            uva = ifft(uafft);
            uvb = ifft(ubfft);
          
            if ((sqrt(norm(uva-u1a,2).^2+norm(uvb-u1b,2).^2) / ...
                 sqrt(norm(u1a,2).^2+norm(u1b,2).^2)) < tol)
            % tolerances met, break loop
                u1a = uva;
                u1b = uvb;
                break;
            else
            % tolerances not met, repeat loop
                u1a = uva;
                u1b = uvb;
            end
        end %end convergence iteration
        if (ii == maxiter)
          warning(sprintf('Failed to converge to %f in %d iterations',...
                          tol,maxiter));
        end
        u0a = u1a;
        u0b = u1b;
    end %end step iteration
  % Convert back from elliptical basis to linear basis:
    u1x = (cos(psi)*cos(chi) + 1j*sin(psi)*sin(chi))*u1a + ...
         (-sin(psi)*cos(chi) - 1j*cos(psi)*sin(chi))*u1b;
    u1y = (sin(psi)*cos(chi) - 1j*cos(psi)*sin(chi))*u1a + ...
          (cos(psi)*cos(chi) - 1j*sin(psi)*sin(chi))*u1b;
% u1x = u1a;
% u1y = u1b;
else
    error('Invalid method specified: %s\n', method);
end
