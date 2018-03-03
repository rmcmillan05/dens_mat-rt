global step
global ck_QS
global E_CS
global E_QS
global E0
global start_delta
global start_delta_id
global V_QS
global ck_CS
global V_CS
global nk_CS
global nk_QS
global step_id
global omegak_CS
global omegak_QS
global gammak_CS
global gammak_QS

out_folder = './OUTPUT/';
out_name = 'test';

max_T = 80; % In femtoseconds
step = 1*1e-3; % In femtoseconds
R_angs = 6.77; % Separation distance in angstrom.
E0=0.0001; % Field intensity in W/cm2

% QS_file is the fitted parameters from the polarizability of the "primary system"
% CS_file is the fitted parameters from the polarizability of the "secondary system"
QS_file = '../data/alpha.dat.fit_params.mat';
CS_file = '../data/alpha.dat.fit_params.mat';

%% Parameters

ha2ev = 27.2113834;
fs2au = 41.3413745758;
angs2bohr = 1.889725989;

%% Get fitting parameters from files
%QS
load(QS_file)

ck_QS = c_k;
nk_QS = s;
omegak_QS = omega_k / ha2ev;
gammak_QS = gamma_k / ha2ev;
ck_QS = ck_QS / ha2ev;

%CS
load(CS_file)

ck_CS = c_k;
nk_CS = s;
omegak_CS = omega_k / ha2ev;
gammak_CS = gamma_k / ha2ev;
ck_CS = ck_CS/ha2ev;

clear('s')
clear('thetak')
clear('omegak')
clear('gammak')

%omegak_CS = omegak_CS + 0.45/ha2ev; % GW corrections
%omegak_QS = omegak_QS + 0.45/ha2ev;

%% Simulation

N = round(max_T/step);

R = R_angs * angs2bohr;

s_alpha = -1; % for X or Y alignment
%s_alpha = 2; % for Z alignment
%s_alpha = -1 * exp(-6.77/3.16); % I've added the exponential factor for the 2D materials.
V_QS = 1;
V_CS = 1;

every_percent = 5;
checkp = N*every_percent/100;
checkt = 0.025*fs2au;

step = step*fs2au;

checkt = round(checkt/step);
start_delta = checkt*step;
start_delta_id = checkt + 1;

t = 0;
t2 = 0;
s_CS = zeros(nk_CS,1);
s_QS = zeros(nk_QS,1);

T = zeros(N,1);
P_QS = zeros(N,1);
P_CS = zeros(N,1);
E = zeros(N,1);
E_CS_out = E;
E_QS_out =E;

for i = 2:N+1
    step_id = i;
    E_CS = eext(t) + s_alpha*P_QS(i-1)/R;
    E_QS = eext(t) + s_alpha*P_CS(i-1)/R;
    E_CS_out(i) = E_CS;
    E_QS_out(i) = E_QS;
    [t, s_CS] = rk_step('sk_CS',t,s_CS);
    [t2, s_QS] = rk_step('sk_QS',t2,s_QS);
    E(i) = eext(t);
    T(i) = t/fs2au;
    P_CS(i) = sum(ck_CS.*real(s_CS));
    P_QS(i) = sum(ck_QS.*real(s_QS));
    if mod(i,checkp) == 0
        fprintf('%.2f %s\n',i/checkp*every_percent,'%')
    end
end

X(:,1) = T(1:checkt:N+1);
X(:,2) = E(1:checkt:N+1)/checkt; %We multiply by this scale to give correct FFT
out_file = [out_folder,out_name,'.dat'];
X(:,3) = P_CS(1:checkt:N+1);
X(:,4) = P_QS(1:checkt:N+1);
X(:,5) = X(:,3) + X(:,4);
X(:,6) = E_CS_out(1:checkt:N+1)/checkt;
X(:,7) = E_QS_out(1:checkt:N+1)/checkt;

if ~isdir(out_folder)
    mkdir(out_folder)
end

save(out_file,'X','-ascii')

fprintf('\n%s%s\n%s\n\n','Output saved to ',out_file,'(Time in column 1, Field in column 2, Polarization (CS) in column 3, Polarization (QS) in column 4, Total Polarization in column 5)')
