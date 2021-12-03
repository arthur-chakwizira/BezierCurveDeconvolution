function [Rut,Psvd,OI] = oSVD_residue_function(C,AIF,time,Psvd,oscFact)

num_samples = max(size(C)); %@@@70


found = 0;
while found == 0
    
    OI = 0; %@@@what is OI?
    
    % run cSVD
    [Rut] = cSVD(C,AIF,time,Psvd); %@@@ residue function calculated using cSVD
    
    %calc Osclillation Index (OI) %@@@Aha, that's OI
    
    %@@@oSVD: Oscillation-index cSVD (oSVD), an iterative method repeating the cSVD deconvolution process
%   until the oscillation in the residue function is below a threshold
%   defined by the oscillation index (OI)
%@@@- cSVD: block-circulant SVD (cSVD), extension of sSVD, which has been shown to be less sensitive to tracer
%delays
%@@@sSVD: standard truncated SVD (sSVD), as originally proposed by Ostergaard et al

    R_scaled = Rut;%./max(Rut);
    
    for k = 3:num_samples %@@@3 to 70
        OI = OI + abs(R_scaled(k)-2*R_scaled(k-1)+R_scaled(k-2)); %@@@ some calculation of OI, will look at it later
    end
    OI = (1/num_samples)*(1/max(Rut))*OI; %@@@ some calculation of OI

    % if Oscillation Index is less than wished-for: exit loop. else increase Psvd and continue
    if OI > oscFact %@@@ if the oscillation index is above oscFact, 
        Psvd = Psvd + 0.02; %increase Psvd   %@@@this means zeroing more elements in the Wc matrix in cSVD.
    else
        found = 1; %exit loop
%         figure(2); hold on
%         plot(time,6000.*Rut,'r'); axis([0 4 -20 80]); title(['Iteration number ' num2str(Psvd*100-19)])
    end

    
    
end