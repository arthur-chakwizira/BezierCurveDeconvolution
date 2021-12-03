function [oef_BzD, cmro2_BzD, oef_svd, cmro2_svd] = oef_cmro2_gui(handles)
%       This function accepts the handles structure and returns OEF and
%       CMRO2, for BzD and SVD (using SVD here is not recommended)
%       Author:
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%display message
% wrap_text = 'Calculating OEF and/or CMRO2...';
% set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
%get needed data__________________________________________________________
mask = handles.mask;
fitd_omega = handles.fitd_omega;
fitd_cbf = handles.fitd_cbf;
fitd_cbf_svd = handles.fitd_cbf_svd;
img_size = handles.img_size;
slice_range = handles.slice_range;
save_results = handles.save_results;
target_folder = handles.target_folder;
fitd_r_svd = handles.fitd_r_svd;
tr = handles.tr;
t = 0:tr:(img_size(4)-1)*tr;
%From Mouridsen 2014: Reliable estimation of capillary transit time
%distributions using DSC-MRI
%dC/dx = -k*tau(a_h*P_50(C/(B-C))^(1/h) - a_h*C_t);
%C is the oxygen concentration and the oxygen extraction fraction for a
%given transit time is given by Q = 1-C(0)/C(1)
%Maximum oxygen extraction is then given by
%OEFmax = integral(0, infinity) h(tau)*Q(tau)dtau
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
    %@@@ integrates the system of differential equations C' = f(C,x) from
    %tspan(0) to tspan(1), with initial conditions given by c_0
    Q(i) = 1 - C(2)/C(1);
end

oef_BzD   = zeros(img_size(1:3));
cmro2_BzD = zeros(img_size(1:3));
oef_svd   = zeros(img_size(1:3));
cmro2_svd = zeros(img_size(1:3));

for x = slice_range(1):slice_range(2)
    for y = slice_range(3):slice_range(4)
        for z = slice_range(5):slice_range(6)
            if mask(x,y,z)
                if handles.BzD
                    r = bezier_residue_function(fitd_omega(x,y,z,:), taus);
                    h = -diff(r)/dtau;
                    h(end+1) = 0;
                    oef_BzD(x,y,z) = trapz(taus,h.*Q);
                    cmro2_BzD(x,y,z) = c_a * oef_BzD(x,y,z) * fitd_cbf(x,y,z);
                end
                if handles.do_SVD
                    try
                        r = interp1(t, fitd_r_svd(x,y,z,:), taus);
                    catch
                        r = zeros(size(taus));
                    end
                    h = -diff(r)/dtau;
                    h(end+1) = 0;
                    oef_svd(x,y,z) = trapz(taus,h.*Q);
                    cmro2_svd(x,y,z) = c_a * oef_svd(x,y,z) * fitd_cbf_svd(x,y,z);
                end
            end
        end
    end
end

%save results
if save_results
    try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
    save_data.this_folder = target_folder; save_data.this_format = handles.save_format;
    if handles.BzD; save_data.this_name = 'oef_BzD'; save_data.data_to_save = oef_BzD; save_this_file(save_data); end
    if handles.do_SVD; save_data.this_name = 'oef_SVD'; save_data.data_to_save = oef_svd; save_this_file(save_data); end
    if handles.BzD; save_data.this_name = 'cmro2_BzD'; save_data.data_to_save = cmro2_BzD; save_this_file(save_data); end
    if handles.do_SVD; save_data.this_name = 'cmro2_SVD'; save_data.data_to_save = cmro2_svd; save_this_file(save_data); end
end


end