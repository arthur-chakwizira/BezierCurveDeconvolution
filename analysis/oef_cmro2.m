function [oef_BzD, cmro2_BzD] = oef_cmro2(data,  config)

cprintf('*blue', 'CALCULATING OEF AND CMRO2 FOR BEZIER CURVE DECONVOLUTION \n')
cprintf('*blue', '\n')

mask = data.mask;
fitd_omega = data.fitd_omega;
fitd_cbf = data.fitd_cbf;
img_size = config.img_size;
slice_range = config.slice_range;
save_results = config.save_results;
target_folder = config.target_folder;

%From Mouridsen 2014: Reliable estimation of capillary transit time
%distributions using DSC-MRI
%dC/dx = -k*tau(a_h*P_50(C/(B-C))^(1/h) - a_h*C_t);
%C is the oxygen conventration and the oxygen extraction fraction for a
%given transit time is given by Q = 1-C(0)/C(1)
%Maximum oxygen extraction is then given by
%OEFmax = integral(0, ?nfinity) h(tau)*Q(tau)dtau
%where h(tau) is the distribution of transit tiems (h(t = )-dR/dt)

k       = 1000;    	% [1/s] permeabiility of capilary wall to O2
%tau is the capillary transit time

% k       = 118;      % [1/s]

B       = 0.1943;   % ml/ml
c_a     = 0.95*B;   % ml/ml
c_0     = c_a;

x   = [0 1];
taus    = (0:0.01:100)';
dtau    = taus(2)-taus(1);

options = odeset('AbsTol',1e-9); %% set solver options
Q = zeros(length(taus),1);

for i = 1:length(taus)
    
    tau = taus(i);
    
    odefun = @(x, C)dc_dx(x,C, tau, k);
    [x,C] = ode45(odefun, x,c_0,options); %% solve equations
    %integrates the system of differential equations C' = f(C,x) from
    %tspan(0) to tspan(1), with initial conditions given by c_0
    Q(i) = 1 - C(2)/C(1);
end

oef_BzD   = zeros(img_size(1:3)); %
cmro2_BzD = zeros(img_size(1:3));

% loop over all voxels
parfor_progress(img_size(1));
for x = slice_range(1):slice_range(2)%1:img_size(1)
    for y = slice_range(3):slice_range(4)%1:img_size(2)
        for z = slice_range(5):slice_range(6)%img_size_3
            if mask(x,y,z)
                r = bezier_residue_function(fitd_omega(x,y,z,:), taus);
                h = -diff(r)/dtau; %distribution of transit times
                h(end+1) = 0;
                oef_BzD(x,y,z) = trapz(taus,h.*Q);
                cmro2_BzD(x,y,z) = c_a * oef_BzD(x,y,z) * fitd_cbf(x,y,z);
            end
        end
    end
    parfor_progress;
end
parfor_progress(0);


if save_results
    niftiwrite(oef_BzD, strcat(target_folder, '/oef_BzD.nii'))
    niftiwrite(cmro2_BzD, strcat(target_folder, '/cmro2_BzD.nii'))
end


cprintf('magenta','Done. \n')
cprintf('blue','_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n')

end