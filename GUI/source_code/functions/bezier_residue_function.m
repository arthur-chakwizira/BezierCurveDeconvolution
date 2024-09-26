function r = bezier_residue_function(omega, t)
%       This function accepts Bezier control points (omega) and a time
%       vector (t), and returns the corresponding Bezier residue function r
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tb = 0:0.01:1; %tau vector for evaluation of Bezier basis polynomials

% if order == 5 %cubic curves
    x = [0 omega(1) omega(3) omega(5)];
    y = [1 omega(2) omega(4)  0    ];
    bx = (1-tb).^3.*x(1) + 3*(1-tb).^2.*tb.*x(2) + 3*(1-tb).*tb.^2.*x(3) + tb.^3.*x(4);
    by = (1-tb).^3.*y(1) + 3*(1-tb).^2.*tb.*y(2) + 3*(1-tb).*tb.^2.*y(3) + tb.^3.*y(4);
% end

% let residue function be zero at the end of the given time vector
bx(end) = t(end);
by(end) = 0;
% resample onto t-grid

try
    r = interp1(bx,by,t,'pchip'); %attempt to evaluate residue function by interpolation
 catch 
     r = t.*0+1; %return ones if that fails
 end
end