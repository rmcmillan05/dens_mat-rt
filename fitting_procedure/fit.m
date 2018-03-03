global s

% File containing the frequency-dependent polarisability, \alpha(\omega), 
% to be fitted
% File format:
%   Column 1: \omega (eV) (Must start from 0eV)
%   Column 2: Imag[\alpha(\omega)]
%   Column 3: Real[\alpha(\omega)]
%
infile = '../data/alpha.dat';

% Number of fitting functions
%
s = 12;

% Limits on parameters for fitting algorithm (units in eV)
%
min_omega_k = -3;
max_omega_k = 50;
min_gamma_k = 0;
max_gamma_k = 10;
min_c_k     = -10;
max_c_k     = Inf;

% Range over which to fit (units in eV)
% If fit_range = [0,0], use entire range available
fit_range = [0,10];

% Required number of points to use in fit. If fit_points = 0 or is less
% than the data available, use all available points
%
fit_points = 0;

% Output file containing the fitted parameters saved as 
% <input_file>.fit_params
%
[pathstr,name,ext] = fileparts(infile);
param_outfile = [infile,'.fit_params'];
% Output file also saved as a MATLAB binary file for easy import back to 
% MATLAB
%
param_outfile_mat = [param_outfile,'.mat'];

% Importing input file
%
fprintf('Reading polarizability from file ''%s''...\n', infile)
fid = fopen(infile);
    A = textscan(fid,'%f%f%f', 'CommentStyle','#');
    imag_col=2;
    real_col=3;
fclose(fid);
omega_orig = A{1};
alpha_orig = (A{real_col} + 1i*A{imag_col});

if fit_range(1)==fit_range(2)
    fit_range = [min(omega_orig), max(omega_orig)];
end

% Selecting data from required fit range
%
fit_range_id = omega_orig >= fit_range(1) & omega_orig <= fit_range(2);
omega_in = omega_orig(fit_range_id);
alpha_in = alpha_orig(fit_range_id);

if fit_points ~= 0
    N = size(omega_in,1);
    fit_points = min(N,fit_points);
    new_range = 1:floor(N/fit_points):N;
    omega_in = omega_in(new_range);
    alpha_in = alpha_in(new_range);
end

N = size(omega_in,1);

% Fitting algorithm requires a minimum value of data points
%
fprintf('%s %i.\n','Number of points fitted = ',N)
if 3*s > 2*N
    error(['3s > 2*N so fit not possible. Either decrease s or ', ...
           'increase the resolution of the input.'])
end
    
% We assume the input file has zomega starting at 0. To fit, we want the
% reflection about the yaxis which should obey the properties of the
% polarisability:
%   Imag[ \alpha(\omega<0) ] = -Imag[ \alpha(\omega>0) ]
%   Real[ \alpha(\omega<0) ] =  Real[ \alpha(\omega>0) ]

%omega = zeros(2*N,1);
%alpha = zeros(2*N,1);

omega(N+1:2*N) = omega_in;      % \omega > 0
omega(1:N) = -omega_in(N:-1:1); % \omega < 0

alpha(N+1:2*N) = alpha_in;
alpha(1:N) = conj(alpha_in(N:-1:1));

% Set up for the least squares procedure
%
v0 = zeros(3*s, 1);
v0(1:s) = [1:s] * max(omega)/(s+1);
v0(s+1:2*s) = 0.025*max(omega)/s;
v0(2*s+1:3*s) = 1.0;

% Lower bounds
%
lb = zeros(3*s, 1);
lb(1:s) = min_omega_k;
lb(s+1:2*s) = min_gamma_k;
lb(2*s+1:3*s) = min_c_k;

% Upper bounds
%
ub = zeros(3*s, 1);
ub(1:s) = max_omega_k;
ub(s+1:2*s) = max_gamma_k;
ub(2*s+1:3*s) = max_c_k;

% We use the trust-region-reflective least squares algorithm
%
opts = optimset('Algorithm', 'trust-region-reflective','Display','Off');

% We use MATLAB's lsqcurvefit to perform a least-squares fit of the 
% input \alpha(\omega) to the fitting functions to obtain the parameters.
% We fit the real and imaginary parts separately by means of the
% function realimag.m as the chosen algorithm is not designed for
% complex functions.
[fitted_parameters,resnorm] = lsqcurvefit(@fitting_function_realimag, v0, omega, ...
                                          realimag(alpha), lb, ub, opts);
omega_k = fitted_parameters(1:s);
gamma_k = fitted_parameters(s+1:2*s);
c_k     = fitted_parameters(2*s+1:3*s);

% This is the approximate \alpha(\omega) from the found parameters
%
alpha_approx = fitting_function(fitted_parameters, omega_orig);

% Estimated percentage error:
%
pe = 100*max(abs(abs(alpha_approx(fit_range_id)) - ...
             abs(alpha_orig(fit_range_id))) / ...
             max(abs(alpha_orig(fit_range_id))));
fprintf('Fitted with estimated percentage error of %.2f%%.\n', pe)  

% Plotting
%
pl = plot(omega_orig, real(alpha_orig), omega_orig, imag(alpha_orig), omega_orig, real(alpha_approx),'-k', omega_orig, imag(alpha_approx),'-k'); 
set(pl(1:2),'linewidth',2)
legend('Real[\alpha]','Imag[\alpha]','Least-Squares Fit')
xlabel('Energy (eV)')
title(['N = ',num2str(s,'%d'),', Fit Range = [',num2str(fit_range(1)),',',num2str(fit_range(2)),']'])

% Saving the paramters to the output files
%
fid = fopen(param_outfile,'w');
fprintf(fid,'%s\n','#');
fprintf(fid,'%s%s\n','# Input file: ',[name, ext]);
fprintf(fid,'%s%f%s%f%s\n','# Fitted over the range: ', ...
                            fit_range(1),':',fit_range(2),' eV');
fprintf(fid,'%s%f%s\n','# Estimated percentage error: ',pe,' %');
fprintf(fid,'%s','#');
fprintf(fid,'\n%s\n', 'n_k');
fprintf(fid,'%d\n', s);
fprintf(fid,'\n%s\n','omega_k');
fprintf(fid,'%f\n',omega_k);
fprintf(fid,'\n%s\n','gamma_k');
fprintf(fid,'%f\n',gamma_k);
fprintf(fid,'\n%s\n','c_k');
fprintf(fid,'%f\n',c_k);
fclose(fid);

fprintf('Fitting parameters saved to file ''%s''.\n', param_outfile)

save(param_outfile_mat,'gamma_k','c_k','omega_k','s')
