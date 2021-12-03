function r = bezier_residue_function(omega, t)
%       This function accepts Bezier control points (omega) and a time
%       vector (t), and returns the corresponding Bezier residue function r
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tb = 0:0.01:1; %tau vector for evaluation of Bezier basis polynomials

order = length(omega);  %how many control points have been supplied

if order == 9 %penta curves
    x = [0 omega(1) omega(3) omega(5) omega(7) omega(9)];
    y = [1 omega(2) omega(4) omega(6)  omega(8)  0     ];
    bx = (1-tb).^5.*x(1) + 5*(1-tb).^4.*tb.*x(2) + 10*(1-tb).^3.*tb.^2.*x(3) + 10*(1-tb).^2.*tb.^3.*x(4) + 5*(1-tb).*tb.^4.*x(5) + tb.^5.*x(6);
    by = (1-tb).^5.*y(1) + 5*(1-tb).^4.*tb.*y(2) + 10*(1-tb).^3.*tb.^2.*y(3) + 10*(1-tb).^2.*tb.^3.*y(4) + 5*(1-tb).*tb.^4.*y(5) + tb.^5.*y(6);
end    

if order == 7   %quartic curves
    x = [0 omega(1) omega(3) omega(5) omega(7)];
    y = [1 omega(2) omega(4) omega(6)   0     ];
    bx = (1-tb).^4.*x(1) + 4*(1-tb).^3.*tb.*x(2) + 6*(1-tb).^2.*tb.^2.*x(3) + 4*(1-tb).*tb.^3.*x(4) + tb.^4.*x(5);
    by = (1-tb).^4.*y(1) + 4*(1-tb).^3.*tb.*y(2) + 6*(1-tb).^2.*tb.^2.*y(3) + 4*(1-tb).*tb.^3.*y(4) + tb.^4.*y(5);
end
if order == 5 %cubic curves
    x = [0 omega(1) omega(3) omega(5)];
    y = [1 omega(2) omega(4)  0    ];
    bx = (1-tb).^3.*x(1) + 3*(1-tb).^2.*tb.*x(2) + 3*(1-tb).*tb.^2.*x(3) + tb.^3.*x(4);
    by = (1-tb).^3.*y(1) + 3*(1-tb).^2.*tb.*y(2) + 3*(1-tb).*tb.^2.*y(3) + tb.^3.*y(4);
end

if order == 3 %quadratic curves
    x = [0 omega(1) omega(3)];
    y = [1 omega(2) 0];
    bx = (1-tb).^2*x(1) + 2*(1-tb).*tb.*x(2) + tb.^2.*x(2);
    by = (1-tb).^2*y(1) + 2*(1-tb).*tb.*y(2) + tb.^2.*y(2);    
end

if order == 4 %cubic curves with ycoordinate of 3rd point fixed
    x = [0 omega(1) omega(3) omega(4)];
    y = [1 omega(2) 0 0];
    bx = (1-tb).^3.*x(1) + 3*(1-tb).^2.*tb.*x(2) + 3*(1-tb).*tb.^2.*x(3) + tb.^3.*x(4);
    by = (1-tb).^3.*y(1) + 3*(1-tb).^2.*tb.*y(2) + 3*(1-tb).*tb.^2.*y(3) + tb.^3.*y(4);
end

%with r(0) free and no dispersion modelling
if order == 6
    x = [0 omega(1) omega(3) omega(5)];
    y = [1 omega(2)  omega(4) omega(6)];
    bx = (1-tb).^3.*x(1) + 3*(1-tb).^2.*tb.*x(2) + 3*(1-tb).*tb.^2.*x(3) + tb.^3.*x(4);
    by = (1-tb).^3.*y(1) + 3*(1-tb).^2.*tb.*y(2) + 3*(1-tb).*tb.^2.*y(3) + tb.^3.*y(4);
end

% let residue function be zero at the end of the given time vector
bx(end+1) = t(end);
by(end+1) = 0;
% resample onto t-grid

try
    r = interp1(bx,by,t,'pchip'); %attempt to evaluate residue function by interpolation
 catch 
     r = t.*0+1; %return ones if that fails
 end
end