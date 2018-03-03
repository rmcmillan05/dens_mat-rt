function s = sk_CS(t,sk)
global gammak_CS
global omegak_CS
global E_CS
global V_CS

s = -(gammak_CS + 1i*omegak_CS).*sk + V_CS*1i*E_CS;