function [Rut] = cSVD_residue_function(C,AIF,time,Psvd)
%        This function accepts a tissue concentration curve, AIF, time
%        vector and Psvd. It returns a residue function obtained using cSVD
%        deconvolution.
%        Contact: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_points = max(size(C));
sizeC = size(C) ;
if sizeC(1)<sizeC(2) 
    C = C'; 
end
sizeAIF = size(AIF); 
if sizeAIF(1)<sizeAIF(2)
    AIF = AIF';
end


smookernel = (1/6).*[1 4 1];

AIFsmoo = conv(AIF,smookernel,'same'); 

C = [C;zeros(num_points,1)];
AIFsmoo = [AIFsmoo;zeros(num_points,1)]; 
deltaT = time(2)-time(1); 
A = tril(toeplitz(AIFsmoo))*deltaT; 

D = zeros(2*num_points,2*num_points); 

for i = 1:2*num_points
    for j = 1:2*num_points
        if j<=i
            D(i,j) = A(i,j);
        else
            D(i,j) = A(2*num_points+i-j,1); 
           
        end
    end
end

[Uc,Sc,Vc] = svd(D);

Wc = diag(diag(1./Sc));

Wc(Sc < Psvd*max(Sc(:))) = 0; 

R = Vc * Wc * Uc' * C; 

Rut = R(1:num_points); 

end
