function f_tmp = fitting_function(parameters, omega)                                                
    
global s
                                                                                
omega_k = parameters(1:s);                                                       
gamma_k = parameters(s+1:2*s);                                                        
c_k     = parameters(2*s+1:3*s);
                                                                                
f_tmp = 0.0;                                                           
for k = 1:s 
    f_tmp = f_tmp + c_k(k)*(1.0./(omega_k(k) - omega - gamma_k(k)*1j) + ...       
                            1.0./(omega_k(k) + omega + gamma_k(k)*1j));       
end