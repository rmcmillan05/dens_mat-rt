function y_out = realimag(y_in)

% This function decomposes a complex variable into a vector containing
% alternations of the real and imaginary parts

y_out = zeros(2*length(y_in),1);
jj = 0;
for ii = 1:length(y_in)
    jj = jj + 1;
    y_out(jj) = real(y_in(ii));
    jj = jj + 1;
    y_out(jj) = imag(y_in(ii));
end