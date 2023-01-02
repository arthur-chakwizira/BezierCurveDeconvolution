function dy_dt = dc_dx(x,C,tau,k)
%       This function is used by the function that calculates OEF and CMRO2
%       'oef_cmro2_gui'. For information and help, contact
%               Arthur Chakwizira
%               arthurchakwizira@gmail.com
%               Medical Radiation Physics, Lund University, Sweden
aH      = 3.1E-5;   % [1/(mm Hg)]
P50     = 26;       % [mm Hg]
B       = 0.1943;   % ml/ml
h       = 2.8;      % [ ]
Pt      = 25;       % [mm Hg]

dy_dt = -k*tau*(aH*P50*(C/(B-C))^(1/h)-aH*Pt);