function R = sSVD_residue_function(C, AIF, time, Psvd)
%        This function accepts a tissue concentration curve, AIF, time
%        vector and Psvd. It returns a residue function obtained using 
%        standard truncated SVD (sSVD) deconvolution
%        Contact: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sizeC = size(C) ;
if sizeC(1)<sizeC(2) 
    C = C'; 
end
sizeAIF = size(AIF); 
if sizeAIF(1)<sizeAIF(2)
    AIF = AIF'; 
end

deltaT = time(2)-time(1);
smookernel = (1/6).*[1 4 1];
AIFsmoo = conv(AIF,smookernel,'same'); % smooth the aif

A = tril(toeplitz(AIFsmoo))*deltaT;
[Uc,Sc,Vc] = svd(A);
Wc = diag(diag(1./Sc)); 
Wc(Sc < Psvd*max(Sc(:))) = 0;

R = Vc * Wc * Uc' * C; 
end