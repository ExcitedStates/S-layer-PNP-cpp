% wriiten by Po-Nan Li @ Stanford 
function [C, C0, c_convs, c_detector, Jx, Jy, Jz, Ex, Ey, Ez, masks] = pnp_car(varargin)

p = inputParser;
p.addParameter('tag', 'cartesian');
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
p.addParameter('slayer', 1);
p.addParameter('save_file', 1);
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
slayer = p.Results.slayer;
save_file = p.Results.slayer;


dx = 1e-10; % Ang

radius = 11e-9; % nm
height = 40e-9; % nm
qpie = dq; % deg
thick_slayer = 4.5e-9;
z_slayer = 20e-9;
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
n_x = floor( radius/dx + 1);
n_y = floor( radius/dx + 1);
n_z = round( height/dx );
n_i = length(charges);
n_t = round( sim_time/dt );


%% PNP-ANA

% ey

sig = ey_sig;
[C_pls, C_min, Ey, ~] = pnp_ana(c_salt, sig, dx, n_z);

C_nh4 = (c_nh4/c_salt) * C_pls;
C_na = (c_salt-c_nh4) / c_salt * C_pls;
C_cl = C_min;

C = c_salt * 1000 * ones( n_z, n_x, n_y, n_i );
C(:,:,:,1) = 1000 * c_nh4;
Jx = zeros( n_z, n_x, n_y, n_i );
Jy = zeros( n_z, n_x, n_y, n_i );
Jz = zeros( n_z+1, n_x, n_y, n_i );
Ex = zeros( n_z, n_x, n_y );
Ey = zeros( n_z, n_x, n_y );
Ez = zeros( n_z+1, n_x, n_y );
c_convs = zeros(1,n_t);
c_detector = zeros(n_i, n_t);


nx_protein = round( r_pore / dx ) + 1;
nz_stop = round( (z_slayer - thick_slayer/2)/dx );
nz_sbtm = round( (z_slayer + thick_slayer/2)/dx );

r = 0:(n_x-1);
[X, Y] = meshgrid(r, r);
R = sqrt(X.^2 + Y.^2);
R = permute(R, [3 1 2]);
for z = 1:n_z
    for k = 1:n_i
        C(z,:,:,k) = C(z,:,:,k) .* (R<n_x);
    end
end
if slayer == 1
    for z = nz_stop:nz_sbtm
        for k = 1:n_i
            C(z,:,:,k) = C(z,:,:,k) .* (R<(nx_protein-1));
        end
    end
end


nz_amo = round(z_amo/dx);

%% J BCs

eps = 1e-6;

Jx_mask = ones(size(Jx));
r = 0:(n_x-1);
[X, Y] = meshgrid(r+0.5, r);
R = sqrt(X.^2 + Y.^2);
R = permute(R, [3 1 2]);
for z = 1:n_z
    for k = 1:n_i
        Jx_mask(z,:,:,k) = Jx_mask(z,:,:,k) .* (R<=(n_x-1+eps));
    end
end

Jy_mask = ones(size(Jy));
r = 0:(n_y-1);
[X, Y] = meshgrid(r, r+0.5);
R = sqrt(X.^2 + Y.^2);
R = permute(R, [3 1 2]);
for z = 1:n_z
    for k = 1:n_i
        Jy_mask(z,:,:,k) = Jy_mask(z,:,:,k) .* (R<=(n_y-1+eps));
    end
end

Jz_mask = ones(size(Jz));

if slayer == 1
    % Jx
    r = 0:(n_x-1);
    [X, Y] = meshgrid(r+0.5, r);
    R = sqrt(X.^2 + Y.^2);
    R = permute(R, [3 1 2]);
    for z = nz_stop:(nz_sbtm)
        for k = 1:n_i
            Jx_mask(z,:,:,k) = Jx_mask(z,:,:,k) .* (R<=(nx_protein-2+eps));
        end
    end
    % Jy
    r = 0:(n_x-1);
    [X, Y] = meshgrid(r, r+0.5);
    R = sqrt(X.^2 + Y.^2);
    R = permute(R, [3 1 2]);
    for z = nz_stop:(nz_sbtm)
        for k = 1:n_i
            Jy_mask(z,:,:,k) = Jy_mask(z,:,:,k) .* (R<=(nx_protein-2+eps));
        end
    end
    % Jz
    r = 0:(n_x-1);
    [X, Y] = meshgrid(r, r);
    R = sqrt(X.^2 + Y.^2);
    R = permute(R, [3 1 2]);
    for z = nz_stop:(nz_sbtm+1)
        for k = 1:n_i
            Jz_mask(z,:,:,k) = Jz_mask(z,:,:,k) .* (R<(nx_protein-1));
        end
    end
end

masks = {};
masks.Jx = Jx_mask;
masks.Jy = Jy_mask;
masks.Jz = Jz_mask;

%% diffusive BCs

C_top = C(1,:,:,:);


%% Run

C0 = C;
filename = [tag '_' int0str(kr,4) '_' int0str(1e4*ex_sig,5) '_' int0str(1e4*ey_sig,5) '_' int0str(1e11*r_pore, 4) '_' int0str(1e9*c_nh4, 6)];

disp(filename);

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
    


tic;


for t = 1:n_t
    
    Cprev = C;
    
    % update E
    for k = 1:n_i
        Ex = Ex - (dt*faraday/eps) * charges(k) * Jx(:,:,:,k);
        Ey = Ey - (dt*faraday/eps) * charges(k) * Jy(:,:,:,k);
        Ez = Ez - (dt*faraday/eps) * charges(k) * Jz(:,:,:,k);
    end

    % update J
    for k = 1:n_i
        Jx(:,1:(end-1),:,k) = -d_m(k)/dx * diff(  C(:,:,:,k), 1, 2 )...
            + 1*(faraday/rc/tmp) * d_m(k) * charges(k) ...
            * 0.5 * ( C(:,1:(end-1),:,k) + C(:,2:end,:,k) ) .* Ex(:,1:(end-1),:);
        Jy(:,:,1:(end-1),k) = -d_m(k)/dx * diff(  C(:,:,:,k), 1, 3 )...
            + 1*(faraday/rc/tmp) * d_m(k) * charges(k) ...
            * 0.5 * ( C(:,:,1:(end-1),k) + C(:,:,2:end,k) ) .* Ey(:,:,1:(end-1));
        Jz(1:(end-1),:,:,k) = -d_m(k)/dx * diff( [C_top(:,:,:,k); C(:,:,:,k)], 1, 1 )...
            + 1*(faraday/rc/tmp) * d_m(k) * charges(k) * ...
            0.5 * ( [C_top(:,:,:,k); C(1:(end-1),:,:,k)] + C(1:end,:,:,k) ) .* Ez(1:(end-1),:,:);
    end

    % BCs for J
    Jx = Jx .* Jx_mask;
    Jy = Jy .* Jy_mask;
    Jz = Jz .* Jz_mask;

    % update C
    for k = 1:n_i
        C(:,:,:,k) = C(:,:,:,k) - dt/dx * diff( cat(2, -Jx(:,1,:,k), Jx(:,:,:,k)), 1, 2 );
        C(:,:,:,k) = C(:,:,:,k) - dt/dx * diff( cat(3, -Jy(:,:,1,k), Jy(:,:,:,k)), 1, 3 );
        C(:,:,:,k) = C(:,:,:,k) - dt/dx * diff( Jz(:,:,:,k), 1, 1 );
    end
    C( C < 0 ) = 0;

    % reaction
    C(nz_amo, 1, 1, 1) = C(nz_amo, 1, 1, 1) - kr * dt / km / (6e23) / dx^3 * C(nz_amo, 1, 1, 1); 

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
if save_file == 1
    save(filename, ...
    'C', 'J*', 'E*', 'C0', 'sim_time', 'dx', 'dt', 'c_convs', 'c_detector');
end

end 
    


