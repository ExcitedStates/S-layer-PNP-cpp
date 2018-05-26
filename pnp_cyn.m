% wriiten by Po-Nan Li @ Stanford for CS279 final project
function [C, C0, c_convs, c_detector, Jr, Jz, Jq, Er, Ez, Eq] = pnp_cyn(varargin)
% sim_time, pnp, movie, kr, ex_sig, ey_sig, r_pore, c_nh4

v3 = 1;
mm = 0;

p = inputParser;
p.addParameter('tag', 'results_pnp_cyn');
p.addParameter('enable_mex', 1);
p.addParameter('n_display', 10000);
p.addParameter('sim_time', 5e-6);
p.addParameter('pnp', 1);
p.addParameter('movie', 0);
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
p.addParameter('z_amo', 30e-9);
p.addParameter('c_limit', 20e-9);
p.parse(varargin{:});

tag = p.Results.tag;
enable_mex = p.Results.enable_mex;
n_display = p.Results.n_display;
sim_time = p.Results.sim_time;
pnp = p.Results.pnp;
movie = p.Results.movie;
kr = p.Results.kr;
ex_sig = p.Results.ex_sig;
ey_sig = p.Results.ey_sig;
r_pore = p.Results.r_pore;
c_nh4 = p.Results.c_nh4;
c_salt = p.Results.c_salt;
cont = p.Results.cont;
dq = p.Results.dq;
dt = p.Results.dt;
z_amo = p.Results.z_amo;
c_limit = p.Results.c_limit;


% dt = 1e-12; % ps
dx = 1e-10; % Ang
% dq = 10; % deg

radius = 11e-9; % nm
height = 40e-9; % nm
qpie = dq; % deg
% z_amo = 30e-9;
% r_pore = 6.5e-10;
thick_slayer = 4.5e-9;
z_slayer = 20e-9;
% kr = 10; % 10 / sec
km = 133e-9 * 1000; 
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
n_r = floor( radius/dx + 1);
n_z = round( height/dx );
n_q = round( qpie/dq );
n_i = length(charges);
n_t = round( sim_time/dt );


%% PNP-ANA

% ey

sig = ey_sig;
[C_pls, C_min, Ey, ~] = pnp_ana(c_salt, sig, dx, n_z);

C_nh4 = (c_nh4/c_salt) * C_pls;
C_na = (c_salt-c_nh4) / c_salt * C_pls;
C_cl = C_min;

C = c_salt * 1000 * ones( n_z, n_r, n_q, n_i );
C(:,:,:,1) = 1000 * c_nh4;
Jr = zeros( n_z, n_r, n_q, n_i );
Jz = zeros( n_z+1, n_r, n_q, n_i );
Jq = zeros( n_z, n_r, n_q+1, n_i );
Er = zeros( n_z, n_r, n_q );
Ez = zeros( n_z+1, n_r, n_q );
Eq = zeros( n_z, n_r, n_q+1 );
c_convs = zeros(1,n_t);
c_detector = zeros(n_i, n_t);


nx_protein = round( r_pore / dx ) + 1;
nz_stop = round( (z_slayer - thick_slayer/2)/dx );
nz_sbtm = round( (z_slayer + thick_slayer/2)/dx );


if pnp > 0 && ey_sig > 0

    C(1:(nz_stop-1),:,:,1) = 1000*repmat( flipud( C_nh4(1:(nz_stop-1)).'), [1 size(C,2) size(C,3) 1] );
    C(1:(nz_stop-1),:,:,2) = 1000*repmat( flipud( C_cl(1:(nz_stop-1)).'), [1 size(C,2) size(C,3) 1] );
    C(1:(nz_stop-1),:,:,3) = 1000*repmat( flipud( C_na(1:(nz_stop-1)).'), [1 size(C,2) size(C,3) 1] );
    C((nz_sbtm+1):end,:,:,1) = 1000*repmat( ( C_nh4(1:(n_z - nz_sbtm)).'), [1 size(C,2) size(C,3) 1] );
    C((nz_sbtm+1):end,:,:,2) = 1000*repmat( ( C_cl(1:(n_z - nz_sbtm)).'), [1 size(C,2) size(C,3) 1] );
    C((nz_sbtm+1):end,:,:,3) = 1000*repmat( ( C_na(1:(n_z - nz_sbtm)).'), [1 size(C,2) size(C,3) 1] );
    C(nz_stop:nz_sbtm, 1:(nx_protein-1), :, 1) = 1000*C_nh4(1);
    C(nz_stop:nz_sbtm, 1:(nx_protein-1), :, 2) = 1000*C_cl(1);
    C(nz_stop:nz_sbtm, 1:(nx_protein-1), :, 3) = 1000*C_na(1);
    
    if v3 == 1
        Ez(1:nz_stop,:,:) = -repmat( flipud(Ey(1:nz_stop).'), [1 size(C,2) size(C,3)] );
        Ez((nz_sbtm+1):end,:,:) = repmat( Ey(1:(n_z-nz_sbtm+1)).', [1 size(C,2) size(C,3)] );
    else
        Ez(1:nz_stop,nx_protein:end,:) = -repmat( flipud(Ey(1:nz_stop).'), [1 size(C,2)-(nx_protein-1) size(C,3)] );
        Ez((nz_sbtm+1):end,nx_protein:end,:) = repmat( Ey(1:(n_z-nz_sbtm+1)).', [1 size(C,2)-(nx_protein-1) size(C,3)] );
    end
end

% ex

sig = ex_sig;
[C_pls, C_min, Ey, ~] = pnp_ana(c_salt, sig, dx, n_z);

C_nh4 = (c_nh4/c_salt) * C_pls;
C_na = (c_salt-c_nh4) / c_salt * C_pls;
C_cl = C_min;


if pnp > 0 && ex_sig > 0
    Er(nz_stop:nz_sbtm, 1:(nx_protein-1), :) = ...
        -repmat( fliplr(Ey(1:(nx_protein-1))), [(nz_sbtm-nz_stop+1) 1 size(C,3)])...
        +repmat( (Ey((nx_protein-1):(2*nx_protein-3))), [(nz_sbtm-nz_stop+1) 1 size(C,3)]);
    C(nz_stop:nz_sbtm, 1:(nx_protein-1), :, 1) = 1000*repmat(fliplr(C_nh4(1:(nx_protein-1))), [(nz_sbtm-nz_stop+1) 1 size(C,3)]);
    C(nz_stop:nz_sbtm, 1:(nx_protein-1), :, 2) = 1000*repmat(fliplr(C_cl(1:(nx_protein-1))), [(nz_sbtm-nz_stop+1) 1 size(C,3)]);
    C(nz_stop:nz_sbtm, 1:(nx_protein-1), :, 3) = 1000*repmat(fliplr(C_na(1:(nx_protein-1))), [(nz_sbtm-nz_stop+1) 1 size(C,3)]);
end

C(nz_stop:nz_sbtm, nx_protein:end, :, :) = 0;


nz_amo = round(z_amo/dx);

if movie == 1
    figure(999);
end

r = 0:(size(C,2)-1);
R = dx * repmat(r, [size(C,1) 1 (size(C,3)+1) 1]);
dqq = pi/180 * dq;

%% BCs

C_top = C(1,:,:,:);


%% ready to run

C0 = C;
filename = [tag '_' int0str(kr,4) '_' int0str(1e4*ex_sig,5) '_' int0str(1e4*ey_sig,5) '_' int0str(1e11*r_pore, 4) '_' int0str(1e9*c_nh4, 6)];

disp(filename);

%% test mat2bin

a = (dt*faraday/eps);
b = faraday/rc/tmp;
c = kr * dt / km / (6e23) / (0.25*dx^3*pi);
        
rst = mat2bin( filename, 0, C, Ez, Er, Eq, Jz, Jr, Jq, ...
                  charges, d_m, dqq, dx, dt, R, ...
                  a, b, c, C0(1,:,:,:), nz_amo, 1, ...
                  [nz_stop, nz_sbtm, nx_protein], n_display);
              
rst

%% a continued run?

if cont && exist([filename '.mat'], 'file') == 2
    sim_time0 = sim_time;
    load(filename);
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
        c = kr * dt / km / (6e23) / (0.25*dx^3*pi);

        [C, Ez, Er, Eq, Jz, Jr, Jq]...
                = mex_pnp3d( C, Ez, Er, Eq, Jz, Jr, Jq, ...
                  charges, d_m, dqq, dx, dt, R, ...
                  a, b, c, ...
                  C0(1,:,:,:), ...
                  nz_amo, 1, [nz_stop, nz_sbtm, nx_protein], ...
                  n_display);
              
        c_dtr(:,tt+tt_offset) = squeeze( C(300,1,1,:) );
        jz_dtr(:,tt+tt_offset) = squeeze( Jz(201,1,1,:) );
              
        t_now = toc;
        t = tt*n_display;
        eta = t_now/t * (n_t-t);
        disp([num2str(round(1000*t/n_t)/10) '% done, ' num2str(round(eta)) 's to go...']);
        
        if mod(t, 500000) == 0
        disp('saving intermediate file');
            save([filename '_t' int0str(t+t_offset, 7)  '.mat'], ...
              'C', 'Ez', 'Er', 'Jz', 'Jr', 'c_dtr', 'jz_dtr', ...
              '-v7.3');
        end  
    end
    
    save(filename, ...
    'C', 'J*', 'E*', 'C0', 'sim_time', 'dx', 'dt', 'c_dtr', 'jz_dtr');
    
else

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
            Jq(:,2:end,:,k) = -d_m(k)/dqq * diff( cat(3, C(:,2:end,end,k), C(:,2:end,:,k), C(:,2:end,1,k)), 1, 3 ) ./ R(:,2:end,:)...
                + 1*(faraday/rc/tmp) * d_m(k) * charges(k) ...
                * 0.5 * ( cat(3, C(:,2:end,end,k), C(:,2:end,:,k)) + cat(3, C(:,2:end,:,k), C(:,2:end,1,k)) ) ...
                .* Eq(:,2:end,:);
        end

        Jr(nz_stop:nz_sbtm, (nx_protein-1), :, :) = 0;
        Jz(nz_stop:(nz_sbtm+1), nx_protein:end, :, :) = 0;

        % update C
        
        % dC^2/dr^2
        C(:,1:end,:,:) = C(:,1:end,:,:) - dt/dx * diff([-Jr(:,1,:,:) Jr(:,1:end,:,:)], 1, 2 );
        % dC/dr
        C(:,1,:,:) =  C(:,1,:,:) -  dt/dx * 2 * Jr(:,1,:,:);
        C(:,2:end,:,:) = C(:,2:end,:,:) - 0.5 * dt * R(:,2:end,2:end).^(-1) .* ( Jr(:,1:(end-1),:,:) + Jr(:,2:end,:,:) );
        % dC^2/dz^2
        C(:,:,:,:) = C(:,:,:,:) - dt/dx * diff( Jz(:,:,:,:), 1, 1 );
        % dC^2/dq^2
        C(:,2:end,:,:) = C(:,2:end,:,:) - dt/dqq * diff( Jq(:,2:end,:,:), 1, 3 ) ./ repmat( R(:,2:end,2:end), [1 1 1 n_i]);

        C( C < 0 ) = 0;

        % reaction
        if mm == 1
            C(nz_amo, 1, :, 1) = C(nz_amo, 1, :, 1) ...
                - abs(kr * dt / (6e23) / (0.25*dx^3*pi) * C(nz_amo, 1, :, 1) ...
                / ( 1 + C(nz_amo, 1, :, 1) / km  ) ); 
        else
            C(nz_amo, 1, :, 1) = C(nz_amo, 1, :, 1) - kr * dt / km / (6e23) / (0.25*dx^3*pi) * C(nz_amo, 1, :, 1); 
        end

        % stats
        c_conv = sum(abs(Cprev(:)-C(:))) / sum(abs(Cprev(:)));
        c_convs(t) = c_conv;
        c_detector(:,t) = 1e6*squeeze( C(nz_amo, 1, 1, :) );
        if movie == 1 && mod(t,10) == 0
            q_slide = 1;
            imagesc( 1e6 * [fliplr(C(:,2:end,q_slide,1)) C(:,:,q_slide,1)] );
            axis image;
            axis off;
            colormap(jet);
            colorbar;
            caxis([0 20]);
            title(['step ' num2str(t) ' / ' num2str(n_t) ', conv = ' num2str(c_conv,2)]);
            pause(0.001);
        end
        dtime = toc;
        eta = dtime / t * (n_t-t);
        if mod(t,n_display) == 0
            disp([num2str(0.1*round(10*eta)) 's to go']);
        end
    end
    save(filename, ...
    'C', 'J*', 'E*', 'C0', 'sim_time', 'dx', 'dt', 'c_convs', 'c_detector');

end 
    


