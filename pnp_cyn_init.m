% wriiten by Po-Nan Li @ Stanford & SLAC
function [C, Jr, Jz, Jq, Er, Ez, Eq, slp] = pnp_cyn_init(varargin)

v3 = 0;

p = inputParser;
p.addParameter('ex_sig', 0);
p.addParameter('ey_sig', 0);
p.addParameter('r_pore', 6.5e-10);
p.addParameter('c_nh4', 10e-9);
p.addParameter('c_salt', 0.14);
p.addParameter('dq', 10);
p.addParameter('dx', 1e-10);
p.parse(varargin{:});

ex_sig = p.Results.ex_sig;
ey_sig = p.Results.ey_sig;
r_pore = p.Results.r_pore;
c_nh4 = p.Results.c_nh4;
c_salt = p.Results.c_salt;
dq = p.Results.dq;
dx = p.Results.dx;
dq = p.Results.dq;


radius = 11e-9; % nm
height = 40e-9; % nm
qpie = dq; % deg
% z_amo = 30e-9;
% r_pore = 6.5e-10;
thick_slayer = 4.5e-9;
z_slayer = 20e-9;

charges = [1 -1 1];


n_r = floor( radius/dx + 1);
n_z = round( height/dx );
n_q = round( qpie/dq );
n_i = length(charges);



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

nx_protein = round( r_pore / dx ) + 1;
nz_stop = round( (z_slayer - thick_slayer/2)/dx );
nz_sbtm = round( (z_slayer + thick_slayer/2)/dx );
slp = [nz_stop nz_sbtm nx_protein];

if ey_sig > 0

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
        Ez(1:nz_stop,nx_protein:end,:)...
            = -repmat( flipud(Ey(1:nz_stop).'), [1 size(C,2)-(nx_protein-1) size(C,3)] );
        Ez((nz_sbtm+1):end,nx_protein:end,:)...
            = repmat( Ey(1:(n_z-nz_sbtm+1)).', [1 size(C,2)-(nx_protein-1) size(C,3)] );
    end
end

% ex

sig = ex_sig;
[C_pls, C_min, Ey, ~] = pnp_ana(c_salt, sig, dx, n_z);

C_nh4 = (c_nh4/c_salt) * C_pls;
C_na = (c_salt-c_nh4) / c_salt * C_pls;
C_cl = C_min;

if ex_sig > 0
    Er(nz_stop:nz_sbtm, 1:(nx_protein-1), :) = ...
        -repmat( fliplr(Ey(1:(nx_protein-1))), [(nz_sbtm-nz_stop+1) 1 size(C,3)])...
        +repmat( (Ey((nx_protein-1):(2*nx_protein-3))), [(nz_sbtm-nz_stop+1) 1 size(C,3)]);
    C(nz_stop:nz_sbtm, 1:(nx_protein-1), :, 1)...
        = 1000*repmat(fliplr(C_nh4(1:(nx_protein-1))), [(nz_sbtm-nz_stop+1) 1 size(C,3)]);
    C(nz_stop:nz_sbtm, 1:(nx_protein-1), :, 2)...
        = 1000*repmat(fliplr(C_cl(1:(nx_protein-1))), [(nz_sbtm-nz_stop+1) 1 size(C,3)]);
    C(nz_stop:nz_sbtm, 1:(nx_protein-1), :, 3)...
        = 1000*repmat(fliplr(C_na(1:(nx_protein-1))), [(nz_sbtm-nz_stop+1) 1 size(C,3)]);
end

C(nz_stop:nz_sbtm, nx_protein:end, :, :) = 0;



