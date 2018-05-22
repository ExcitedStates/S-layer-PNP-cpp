function rst = mat2bin(tag, t, C, Ez, Er, Eq, Jz, Jr, Jq, charges, d_m, ... 
    dqq, dx, dt, R, a, b, c,bulk, nz_amo, nr_amo, slp, n_display)

rst = 1;

%% init fid

filename = [tag '_t' int0str(t, 9) '.pnl'];
fid = fopen(filename, 'wb');
fwrite(fid, 1, 'int'); % input mode
fwrite(fid, t, 'int'); % input mode

%% write shape

[m0, n0, l0, ~] = size(C);
[m1, n1, l1, ~] = size(Jz);
[m2, n2, l2, ~] = size(Jr);
[m3, n3, l3, ~] = size(Jq);
fwrite(fid, [m0 n0 l0 m1 n1 l1 m2 n2 l2 m3 n3 l3], 'int');

%% write C profile

fwrite(fid, C, 'double');

%% write E profiles

fwrite(fid, Ez, 'double');
fwrite(fid, Er, 'double');
fwrite(fid, Eq, 'double');

%% write J profiles

fwrite(fid, Jz, 'double');
fwrite(fid, Jr, 'double');
fwrite(fid, Jq, 'double');

%% write constants

data = [charges, d_m, dqq, dx, dt];
fwrite(fid, data, 'double');

%% write R

fwrite(fid, R, 'double');

%% write coefficients

data = [a, b, c];
fwrite(fid, data, 'double');

%% write bulk

fwrite(fid, bulk, 'double');

%% write misc.

data = [nz_amo, nr_amo, slp, n_display];
fwrite(fid, data, 'double');

%% close file

fclose(fid);
rst = 0;

end
