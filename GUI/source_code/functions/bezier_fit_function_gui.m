function J = bezier_fit_function_gui(p,t,ydata,para, Pnorm, noise,x0)
%       This is the objective function for Bezier curve deconvolution. It
%       returns a scalar, which is the mode of the posterior distribution
%       of the parameters (BzD is implemented in a Bayesian framework 
%       adopting the Maximum Aposteriori Approach).
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
te = para.te;
tmp_baseline_idx = para.baseline_index;
c_aif       = para.c_aif; 
aif_tail_c  = para.aif_tail_c;
kappa       = para.kappa;
with_dispersion = para.with_dispersion;
with_delay = para.with_delay;

p = exp(p.*Pnorm);
dt = t(2)-t(1);
M = length(t);

%constraints
lambda1 = 0; %use this to penalise undesired outcomes

omega = p(1:5); % using order 3 curves
if (omega(2)>1); lambda1 = 1E6; end
if (omega(4)>1); lambda1 = 1E6; end

if omega(1)>omega(5); lambda1 = 1E6; end
if omega(3)>omega(5); lambda1 = 1E6; end

omega_s       = [   8   1   4  1   100  ]; %prior std of control points

f = p(6); %cbf

if f < 0; lambda1 = 1E6; end %reject negative cbf

f_s = 1E6;   %prior std of cbf; cbf prior is non-informative
omega_0       = [     8   0.5   2  0.0   15   ]; %prior mean of ctrl pnts
f_0 = 0.01; %prior mean of cbf
% -------------------------------------------------------------------------
if with_dispersion && ~with_delay
    dkp = p(7:8); %
    if (dkp(1)<=0); dkp(1) = 1E-3; end
    if (dkp(2)<0); dkp(2) = 1E-3; end
    p(7:8) = dkp;
    dk_s  = [     2      2 ]; %s p 
    %gamma dispersion kernel
    dks   = exp(dkp(1));
    dkp   = exp(dkp(2)); 
    dk  = ((dks^(1+dks*dkp))/gamma(1+dks*dkp)).*t.^(dks*dkp).*exp(-dks.*t);
    delay_s = [];
    delay_0 = [];
    dk_0 = x0(7:8);
    mean_aif_c  = dt*filter(c_aif,1,dk); %disperse aif       
end
% -------------------------------------------------------------------------

if with_delay && ~with_dispersion
    delay = p(7);
    if (delay < -t(end)); delay = 1E-3; end
    if (delay > t(end)); delay = 1E-3; end
    %update p
     p(7) = delay;
    %priors
    delay_s  = 5;
    dk_s = [];
    dk_0 = [];
    delay_0 = x0(7);
    mean_aif_c = interp1(t,c_aif,t-delay,'pchip','extrap'); %delay aif
    if (delay < 0); id = ceil(-delay/dt); 
        mean_aif_c(end-id:end) = aif_tail_c;
    end
    if (delay > 0); id = ceil(delay/dt); 
        mean_aif_c(1:id) = 0; 
    end
end
% -------------------------------------------------------------------------
if with_dispersion && with_delay
    dkp = p(8:9); %
    if (dkp(1)<=0); dkp(1) = 1E-3; end
    if (dkp(2)<0); dkp(2) = 1E-3; end
    p(8:9) = dkp;
    dk_s  = [     2      2 ]; %s p
    %gamma dispersion kernel
    dks   = exp(dkp(1));
    dkp   = exp(dkp(2));  
    dk  = ((dks^(1+dks*dkp))/gamma(1+dks*dkp)).*t.^(dks*dkp).*exp(-dks.*t);
    delay = p(7);
    if (delay < -t(end)); delay = 1E-3; end
    if (delay > t(end)); delay = 1E-3; end
    %update p
     p(7) = delay;
    %priors
    delay_s  = 5;
    delay_0 = x0(7);
    dk_0 = x0(8:9);
    
    mean_aif_c = interp1(t,c_aif,t-delay,'pchip','extrap'); %delay aif
    if (delay < 0); id = ceil(-delay/dt); 
        mean_aif_c(end-id:end) = aif_tail_c;
    end
    if (delay > 0); id = ceil(delay/dt); 
        mean_aif_c(1:id) = 0; 
    end
    
    mean_aif_c  = dt*filter(mean_aif_c,1,dk); %disperse delayed aif   
end
% -------------------------------------------------------------------------

if ~with_dispersion && ~with_delay
    delay_s  = [];
    dk_s = [];
    delay_0 = [];
    dk_0 = [];
    mean_aif_c = c_aif; %do nothing to the aif
end
% -------------------------------------------------------------------------

mux = [omega_0, f_0, delay_0, dk_0];
sigmax = [ omega_s    f_s    delay_s    dk_s ]; %std of prior distribution

if (~isempty(omega(omega<0))); lambda1 = 1E9; end

% bezier residue function
r = bezier_residue_function(omega, t);

%convolve (delayed/dispersed) aif with residue function and convert to
%signal
y_m = kappa*f*dt.*filter(mean_aif_c,1,r);
tissue_s0 = mean(ydata(tmp_baseline_idx));
y_m = tissue_s0.*exp(-te*y_m);

% now we calculate the negative logarithm of the posterior probability distribution, J
J = M * log(noise) + 0.5*(  sum(((ydata-y_m)./noise).^2) + sum(((p-mux)./sigmax).^2) ) + lambda1;
% if ~isfinite(J); J = 1e9;  end
% J = y_m-y_data; %trying lsqnonlin

end
