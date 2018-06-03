% wriiten by Po-Nan Li @ Stanford for CS279 final project
function [C, C0, c_dtr, jz_dtr, c_detector, Jr, Jz, Jq, Er, Ez, Eq] = pnp_cyn(varargin)
% sim_time, pnp, movie, kr, ex_sig, ey_sig, r_pore, c_nh4

v3 = 0;
mm = 1;


p = inputParser;
p.addParameter('tag', 'test');
p.addParameter('enable_mex', 0);
p.addParameter('n_display', 1);
p.addParameter('sim_time', 5e-6);
p.addParameter('kr', 10);
p.addParameter('ex_sig', 0);
p.addParameter('ey_sig', 0);
p.addParameter('r_pore', 6.5e-10);
p.addParameter('c_nh4', 10e-9);
p.addParameter('c_salt', 0.14);
p.addParameter('sig', 0);
p.addParameter('cont', 0);
p.addParameter('dq', 10);
p.addParameter('dt', 1e-12);
p.addParameter('dx', 1e-10);
p.addParameter('z_amo', 30e-9);
p.addParameter('late_start', 0);
p.addParameter('kd', 133e-9);
p.addParameter('make_pnl', 0);
p.addParameter('save_intermediate', 0);
p.parse(varargin{:});

tag = p.Results.tag;
enable_mex = p.Results.enable_mex;
n_display = p.Results.n_display;
sim_time = p.Results.sim_time;
kr = p.Results.kr;
ex_sig = p.Results.ex_sig;
ey_sig = p.Results.ey_sig;
r_pore = p.Results.r_pore;
c_nh4 = p.Results.c_nh4;
c_salt = p.Results.c_salt;
cont = p.Results.cont;
dq = p.Results.dq;
dt = p.Results.dt;
dx = p.Results.dx;
z_amo = p.Results.z_amo;
late_start = p.Results.late_start;
kd = p.Results.kd;
make_pnl = p.Results.make_pnl;
save_intermediate = p.Results.save_intermediate;

% dt = 1e-12; % ps
% dx = 1e-10; % Ang
% dq = 10; % deg
height = 40e-9; % nm
% z_amo = 30e-9;
% r_pore = 6.5e-10;
% kr = 10; % 10 / sec
km = 133e-9 * 1000; 
kd = kd * 1000;
e_c = 1.602e-19; % elementary charge in C
avo = 6.022e23;
rc = 8.31446; % reaction constant
faraday = avo * e_c;
tmp = 300; % Kelvin
eps =  80 * 8.854e-12; % solvent permittivity in C/V/m

charges = [1 -1 1];
d_nh4 = 1.64e-9; % m^2/s
d_cl  = 2.03e-9; % m^2/s
d_na  = 1.33e-9; % m^2/s
d_nacl = 1.61e-9; % m^2/s

d_m = [d_nh4 d_cl d_na];
n_z = round( height/dx );
n_i = length(charges);
n_t = round( sim_time/dt );

t_after = round(late_start / dt);

%% filename

filename = [tag '_' num2sci(kr) '_' num2sci(ex_sig) '_' ...
    num2sci(ey_sig) '_' num2sci(r_pore) '_' ...
    num2sci(c_nh4) '_' num2sci(1e-3*kd)];
disp(filename);

%% get initial profiles

[C, Jr, Jz, Jq, Er, Ez, Eq, slp] = pnp_cyn_init('ex_sig', ex_sig, ...
    'ey_sig', ey_sig, 'r_pore', r_pore, 'c_nh4', c_nh4, 'c_salt', c_salt, ...
    'dq', dq, 'dx', dx);
if cont && exist([filename '.mat'], 'file') == 2
    sim_time0 = sim_time;
    load(filename, 'C', 'J*', 'E*', 'C0', 'sim_time', 'dx', 'dt', 'c_convs', 'c_detector');
    sim_time = sim_time + sim_time0;
    disp([filename ' loaded']);
    disp([num2str(sim_time-sim_time0) 's ran, ' num2str(sim_time0) 's to go']);
    t_offset = round((sim_time-sim_time0)/dt);
    tt_offset = round(t_offset/n_display);
else
    t_offset = 0;
    tt_offset = 0;
    cont = 0;
end
    
nz_amo = round(z_amo/dx);
r = 0:(size(C,2)-1);
R = dx * repmat(r, [size(C,1) 1 (size(C,3)+1) 1]);
dqq = pi/180 * dq;

%% BCs

C_top = C(1,:,:,:);
C0 = C;

%% coefficients for MEX or C++ solver

a = (dt*faraday/eps);
b = faraday/rc/tmp;
c = kr * dt / kd  / (6e23) / (0.25*dx^3*pi);

if make_pnl
    mat2bin( ['../pnp12/' filename], 0, C, Ez, Er, Eq, Jz, Jr, Jq, ...
        charges, d_m, dqq, dx, dt, R, a, b, c, km, C0(1,:,:,:), nz_amo, 1, ...
        slp, n_display, t_after);
    return;
end
              
%% solver

tic;

if enable_mex
    
    n_tt = ceil(n_t/n_display);
    
    if cont == 1
        try 
            c_dtr = [c_dtr zeros(3, n_tt)];
            jz_dtr = [jz_dtr zeros(3, n_tt)];
        catch 
            c_dtr = zeros(3, ceil(sim_time/dt));
            jz_dtr = zeros(3, ceil(sim_time/dt));
        end
    else
        c_dtr = zeros(3, n_tt);
        jz_dtr = zeros(3, n_tt);
    end
    
 
    for tt = 1:n_tt
    
        a = (dt*faraday/eps);
        b = faraday/rc/tmp;
        c = kr * dt / (6e23) / (0.25*dx^3*pi) / kd;

        [C, Ez, Er, Eq, Jz, Jr, Jq]...
                = mex_pnp3d( C, Ez, Er, Eq, Jz, Jr, Jq, ...
                  charges, d_m, dqq, dx, dt, R, ...
                  a, b, c, km, ...
                  C0(1,:,:,:), ...
                  nz_amo, 1, slp, ...
                  n_display, t_after);
              
        t_after = t_after - n_display;     
        c_dtr(:,tt+tt_offset) = squeeze( C(nz_amo,1,1,:) );
        jz_dtr(:,tt+tt_offset) = squeeze( Jz(round((n_z+1)/2),1,1,:) );
              
        t_now = toc;
        t = tt*n_display;
        eta = t_now/t * (n_t-t);
        disp([num2str(round(1000*t/n_t)/10) '% done, ' num2str(round(eta)) 's to go...']);
        
        if mod(t, save_intermediate) == 0
            disp('saving intermediate file');
                save([filename '_t' int0str(t+t_offset, 7)  '.mat'], ...
                  'C', 'Ez', 'Er', 'Jz', 'Jr', 'c_dtr', 'jz_dtr', '-v7.3');
        end  
    end
    
    save([filename '.mat'], ...
    'C', 'J*', 'E*', 'C0', 'sim_time', 'dx', 'dt', 'c_dtr', 'jz_dtr');
    
else
    nz_stop = slp(1);
    nz_sbtm = slp(2);
    nx_protein = slp(3);
    for t = 1:n_t
        
        Cprev = C;
        
        % update E
        for k = 1:n_i
            Er = Er - (dt*faraday/eps) * charges(k) * Jr(:,:,:,k);
            Ez = Ez - (dt*faraday/eps) * charges(k) * Jz(:,:,:,k);
            Eq = Eq - (dt*faraday/eps) * charges(k) * Jq(:,:,:,k);
        end

        % update J
        for k = 1:n_i
            Jr(:,1:(end-1),:,k) = -d_m(k)/dx * diff(  C(:,:,:,k), 1, 2 )...
                + 1*(faraday/rc/tmp) * d_m(k) * charges(k) ...
                * 0.5 * ( C(:,1:(end-1),:,k) + C(:,2:end,:,k) ) .* Er(:,1:(end-1),:);
            Jz(1:(end-1),:,:,k) = -d_m(k)/dx * diff( [C_top(:,:,:,k); C(:,:,:,k)], 1, 1 )...
                + 1*(faraday/rc/tmp) * d_m(k) * charges(k) * ...
                0.5 * ( [C_top(:,:,:,k); C(1:(end-1),:,:,k)] + C(1:end,:,:,k) ) .* Ez(1:(end-1),:,:);
            Jq(:,2:end,:,k) = ...
                -d_m(k)/dqq * diff( cat(3, C(:,2:end,end,k), C(:,2:end,:,k), C(:,2:end,1,k)), 1, 3 ) ./ R(:,2:end,:)...
                + 1*(faraday/rc/tmp) * d_m(k) * charges(k) ...
                * 0.5 * ( cat(3, C(:,2:end,end,k), C(:,2:end,:,k)) + cat(3, C(:,2:end,:,k), C(:,2:end,1,k)) ) ...
                .* Eq(:,2:end,:);
        end

        Jr(nz_stop:nz_sbtm, (nx_protein-1), :, :) = 0;
        Jz(nz_stop:(nz_sbtm+1), nx_protein:end, :, :) = 0;

        % update C
        
        % dC^2/dr^2
        C(:,1:end,:,:) = ...
            C(:,1:end,:,:) - dt/dx * diff([-Jr(:,1,:,:) Jr(:,1:end,:,:)], 1, 2 );
        % dC/dr
        C(:,1,:,:) =  C(:,1,:,:) -  dt/dx * 2 * Jr(:,1,:,:);
        C(:,2:end,:,:) = ...
            C(:,2:end,:,:) - 0.5 * dt * R(:,2:end,2:end).^(-1) ...
            .* ( Jr(:,1:(end-1),:,:) + Jr(:,2:end,:,:) );
        % dC^2/dz^2
        C(:,:,:,:) = C(:,:,:,:) - dt/dx * diff( Jz(:,:,:,:), 1, 1 );
        % dC^2/dq^2
        C(:,2:end,:,:) = ...
            C(:,2:end,:,:) - dt/dqq * diff( Jq(:,2:end,:,:), 1, 3 ) ...
            ./ repmat( R(:,2:end,2:end), [1 1 1 n_i]);

        C( C < 0 ) = 0;

        % reaction
        if t > t_after
            if mm == 1
                C(nz_amo, 1, :, 1) = C(nz_amo, 1, :, 1) ...
                    - abs(kr/kd * dt / (6e23) / (0.25*dx^3*pi) * C(nz_amo, 1, :, 1) ...
                    / ( 1 + C(nz_amo, 1, :, 1) / km  ) ); 
            else
                C(nz_amo, 1, :, 1) = ...
                    C(nz_amo, 1, :, 1) - kr * dt / km / (6e23) / (0.25*dx^3*pi) * C(nz_amo, 1, :, 1); 
            end
        end

        % stats
        c_conv = sum(abs(Cprev(:)-C(:))) / sum(abs(Cprev(:)));
        c_convs(t) = c_conv;
        c_detector(:,t) = 1e6*squeeze( C(nz_amo, 1, 1, :) );

        dtime = toc;
        eta = dtime / t * (n_t-t);
        if mod(t,n_display) == 0
            disp([num2str(0.1*round(10*eta)) 's to go']);
        end
    end
    save([filename '.mat'], ...
    'C', 'J*', 'E*', 'C0', 'sim_time', 'dx', 'dt', 'c_convs', 'c_detector');

end 
    


