function Rut = oSVD_residue_function(C,AIF,time,Psvd, oscFact)
%        This function accepts a tissue concentration curve, AIF, time
%        vector and Psvd. It returns a residue function obtained using oSVD
%        deconvolution.
%        Contact: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_samples = max(size(C)); 
found = 0;
while found == 0
    OI = 0;
    % run cSVD
    [Rut] = cSVD_residue_function(C,AIF,time,Psvd); 
    %calculate Osclillation Index (OI)
    R_scaled = Rut;%./max(Rut);
    for k = 3:num_samples
        OI = OI + abs(R_scaled(k)-2*R_scaled(k-1)+R_scaled(k-2)); 
    end
    OI = (1/num_samples)*(1/max(Rut))*OI; 
    % if Oscillation Index is less than wished-for: exit loop. else increase Psvd and continue
    if OI > oscFact 
        Psvd = Psvd + 0.05; %the increment here has a significant bearing on execution time 
    else
        found = 1; %exit loop
    end  
end