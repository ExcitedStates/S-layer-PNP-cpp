function [C, Ez, Er, Eq, Jz, Jr, Jq, charges, d_m, ... 
    dqq, dx, dt, R, a, b, c,bulk, nz_amo, nr_amo, slp, n_display] = bin2mat(filename)

n_i = 3;

%% init fid

fid = fopen(filename, 'r');
mode = fread(fid, 1, 'int');
t = fread(fid, 1, 'int');

if mode == 1
    disp('mode: input');
else
    disp('mode: output');
end

disp(['from time step ' int2str(t)]);

%% read shapes

shapes = fread(fid, 12, 'int');
shapes0 = shapes(1:3).';
shapes1 = shapes(4:6).';
shapes2 = shapes(7:9).';
shapes3 = shapes(10:12).';

%% read C profile

C = fread(fid, prod([shapes0 n_i]), 'double');
C = reshape( C, [shapes0 n_i] );

%% read E profiles

Ez = fread(fid, prod(shapes1), 'double');
Ez = reshape( Ez, shapes1 );

Er = fread(fid, prod(shapes2), 'double');
Er = reshape( Er, shapes2 );

Eq = fread(fid, prod(shapes3), 'double');
Eq = reshape( Eq, shapes3 );

%% read J profiles

Jz = fread(fid, prod([shapes1 n_i]), 'double');
Jz = reshape( Jz, [shapes1 n_i] );

Jr = fread(fid, prod([shapes2 n_i]), 'double');
Jr = reshape( Jr, [shapes2 n_i] );

Jq = fread(fid, prod([shapes3 n_i]), 'double');
Jq = reshape( Jq, [shapes3 n_i] );

%% early return

if mode == 0
    return
end

%% read constants

data = fread(fid, 9, 'double');
charges = data(1:3);
d_m = data(4:6);
dqq = data(7);
dx = data(8);
dt = data(9);

%% read R

R = fread(fid, prod([shapes0(1:2) (shapes0(3)+1)]), 'double');
R = reshape( R, [shapes0(1:2) (shapes0(3)+1)]);

%% read coefficients

data = fread(fid, 3, 'double');
a = data(1);
b = data(2);
c = data(3);

%% write bulk

bulk = fread(fid, prod([1 shapes0(2:3) n_i]), 'double');
bulk = reshape( bulk, [1 shapes0(2:3) n_i] );

%% write misc.

data = fread(fid, 6, 'double');
nz_amo = data(1);
nr_amo = data(2);
slp = data(3:5);
n_display = data(6);

%% close file

fclose(fid);


end
