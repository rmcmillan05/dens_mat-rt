function s = sk_QS(t,sk)
global gammak_QS
global omegak_QS
global E_QS
global V_QS

s = -(gammak_QS + 1i*omegak_QS).*sk + V_QS*1i*E_QS;