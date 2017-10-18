function [C, C0, c_dtr, jy_dtr] = pnp_num(sig1, kk, cc1, cc2, pore_size, amo_y, amo_x)


label = '50ms';
enable_sym = 1;
enable_mex = 1;

%% quick setting

switch nargin
    case 0
        sig1 = 1;
        kk = 0;
        cc1 = 10e-9;
        cc2 = 0.14;
        pore_size = 1.3;
        amo_y = 555;
        amo_x = 111;
    case 1
        kk = 0;
        cc1 = 10e-9;
        cc2 = 0.14;
        pore_size = 1.3;
        amo_y = 555;
        amo_x = 111;
    case 2
        cc1 = 10e-9;
        cc2 = 0.14;
        pore_size = 1.3;
        amo_y = 555;
        amo_x = 111;
    case 3
        cc2 = 0.14;
        pore_size = 1.3;
        amo_y = 555;
        amo_x = 111;
    case 4
        pore_size = 1.3;
        amo_y = 555;
        amo_x = 111;
    case 5
        amo_y = 555;
        amo_x = 111;
    case 6
        amo_x = 111;
end

sig2 = sig1;

cc = cc1 * 1e9;

if enable_sym
    disp('symmetry is enabled!');
end

disp(['(' num2str(sig1) 'e per A2, ' num2str(kk) ' rate, '...
     num2str(cc) ' nM, ' num2str(cc2) ' M, ' num2str(pore_size)...
     ' nm, amo_x ' num2str(amo_x) ' A, amo_y ' num2str(amo_y) ' A'  ' )']);
        

%% sim setting

close all;

real_time_movie = 0;
reaction = 0;

n_display = 10000;

%% constants

k_b = 1.38064852e-23; % Boltzmann constant in  m2 kg s-2 K-1
e_c = 1.602e-19; % elementary charge in C
avo = 6.022e23;
r = 8.31446;
faraday = avo * e_c;


%% parameters

% simulation
dt = 1e-12; % s
dx = 0.1e-9; % m

l_t = 50000000 * dt; % 1e4 * dt;
l_x = 22.1e-9;
l_y = 60.1e-9;

% experiment

kr = kk; % 12.8 day^-1
d_nh4 = 1.64e-9; % m^2/s
d_cl = 2.03e-9;
d_na = 1.33e-9;
d_nacl = 1.61e-9; % m^2/s
d_m = [d_nh4 d_cl d_na];
eps =  80 * 8.854e-12; % solvent permittivity in C/V/m
eps_o =  80 * 8.854e-12;
tmp = 300;
z = [1 -1 1]; % valence charge
k_m = 133e-9 * 1000; 

%% membrane


x_s = 1;

y_pc = (l_y - 20e-9/x_s); % center of the protein is 20 nm above the bottom

% l_px = 1.3e-9;
l_px = pore_size *1e-9;
l_py = 4.5e-9;

l_px = l_px / x_s;
l_py = l_py / x_s;

if l_px > l_x
    no_slayer = 1;
else
    no_slayer = 0;
end

if no_slayer
    y_pt = [];
    y_pb = [];
    x_pl = [];
    x_pr = [];
else
    y_pt = round( (y_pc-l_py/2)/dx );
    y_pb = round( (y_pc+l_py/2)/dx );
    x_pl = round( (l_x-l_px)/2/dx );
    x_pr = round( ((l_x-l_px)/2+l_px)/dx ) + 1;
end

% domain 1
% C( y_pt:y_pb, 1:end, : ) = 0;
% domain 2
% C( y_pt:y_pb, 1:end, : ) = 0;
%% init

n_i = length(z);

n_t = ceil(l_t/dt);
n_x = ceil(l_x/dx);
n_y = ceil(l_y/dx);

x_c = round( (n_x+1)/2 );
y_c = round( (n_y+1)/2 );

% analytical

c_salt = cc2;

% [C1_pls, C1_min, Ey10, V1] = pnp_ana(cc*1e-9, sig1, dx, n_y);
% [C2_pls, C2_min, Ey20, V2] = pnp_ana(cc*1e-9, sig2, dx, n_y);
[C3_pls, C3_min, Ey30, V3] = pnp_ana(c_salt, sig1, dx, n_y);
[C4_pls, C4_min, Ey40, V4] = pnp_ana(c_salt, sig2, dx, n_y);

C3_pls = (cc2-cc1)/cc2 * C3_pls;
C4_pls = (cc2-cc1)/cc2 * C4_pls;

C1_pls = (cc*1e-9/c_salt) * C3_pls;
C2_pls = (cc*1e-9/c_salt) * C4_pls;

% conecentration
C = zeros( n_y, n_x, n_i);
if sig1 > 0 && ~no_slayer
    C(1:(y_pt-1),:,1) = 1000*repmat( fliplr(C1_pls(2:(y_pt-0))).', [1 n_x]);
    C(1:(y_pt-1),:,2) = 1000*repmat( fliplr(C3_min(2:(y_pt-0))).', [1 n_x]);
    C(1:(y_pt-1),:,3) = 1000*repmat( fliplr(C3_pls(2:(y_pt-0))).', [1 n_x]);
    grad1 = linspace(1, 0, length(y_pt:y_pb)).';
    grad2 = linspace(0, 1, length(y_pt:y_pb)).';
    % ion 1
    C(y_pt:y_pb,(x_pl+1):(x_pr-1),1) ...
      = 1000*( grad1 * C1_pls(2:(length((x_pl+1):(x_pr-1))+1)) + grad2 * (C2_pls(2:(length((x_pl+1):(x_pr-1))+1))) );
    C(y_pt:y_pb,(x_pl+1):(x_pr-1),1) = 0.5 * ( C(y_pt:y_pb,(x_pl+1):(x_pr-1),1) + fliplr( C(y_pt:y_pb,(x_pl+1):(x_pr-1),1) ) );
    % ion 2
    C(y_pt:y_pb,(x_pl+1):(x_pr-1),2) ...
      = 1000*( grad1 * C3_min(2:(length((x_pl+1):(x_pr-1))+1)) + grad2 * (C4_min(2:(length((x_pl+1):(x_pr-1))+1))) );
    C(y_pt:y_pb,(x_pl+1):(x_pr-1),2) = 0.5 * ( C(y_pt:y_pb,(x_pl+1):(x_pr-1),2) + fliplr( C(y_pt:y_pb,(x_pl+1):(x_pr-1),2) ) );
    % ion 3
    C(y_pt:y_pb,(x_pl+1):(x_pr-1),3) ...
      = 1000*( grad1 * C3_pls(2:(length((x_pl+1):(x_pr-1))+1)) + grad2 * (C4_pls(2:(length((x_pl+1):(x_pr-1))+1))) );
    C(y_pt:y_pb,(x_pl+1):(x_pr-1),3) = 0.5 * ( C(y_pt:y_pb,(x_pl+1):(x_pr-1),3) + fliplr( C(y_pt:y_pb,(x_pl+1):(x_pr-1),3) ) );

    C((y_pb+1):end,:,1) = 1000*repmat( (C2_pls(2:(size(C((y_pb+1):end,:,1),1)+1))).', [1 n_x]);
    C((y_pb+1):end,:,2) = 1000*repmat( (C4_min(2:(size(C((y_pb+1):end,:,2),1)+1))).', [1 n_x]);
    C((y_pb+1):end,:,3) = 1000*repmat( (C4_pls(2:(size(C((y_pb+1):end,:,2),1)+1))).', [1 n_x]);
else
    C(:,:,1) = 1e-6 * cc;
    C(:,:,2) = 1000 * c_salt;
    C(:,:,3) = 1000 * c_salt;
    
    % domain 1
    C( y_pt:y_pb, 1:x_pl, : ) = 0;
    % domain 2
    C( y_pt:y_pb, x_pr:end, : ) = 0;
    
end


% C(:,:,2) = 1000*repmat( fliplr(C_min(2:(end-1))).', [1 n_x]);

% C(C~=0) = 10e-9 * 1000;

c_bulk_pls = C(1,:,1);
c_salt_min = C(1,:,2);
c_salt_pls = C(1,:,3);

% C(2:end,:,1) = c_bulk_pls;
% C(2:end,:,2) = c_bulk_min;

% fluxes
Jx = zeros( n_y, n_x+1, n_i);
Jy = zeros( n_y+1, n_x, n_i);
% E-fields
Ex = zeros( n_y, n_x+1);
Ey = zeros( n_y+1, n_x);
% Ey = -1000 * repmat(fliplr(Ey0).', [1 n_x]);
% Ey(1:(end-1),:) = 0;

% ICs

% C(1, :) = 100;
% C(131,107) = 1e0;




%% generate E-fields profile

rho = zeros(size(C(:,:,1)));

sigma = -e_c /25; % custom

rho( y_pt, 1:x_pl ) = 1;
rho( y_pb, 1:x_pl ) = 1;
rho( y_pt, x_pr:end ) = 1;
rho( y_pb, x_pr:end ) = 1;
rho( y_pt:y_pb, x_pl ) = 1;
rho( y_pt:y_pb, x_pr ) = 1;


% [X1, Y1] = meshgrid((0.5*dx):dx:(l_x+0.5*dx), dx:dx:l_y);
% [X2, Y2] = meshgrid(dx:dx:l_x, (0.5*dx):dx:(l_y+0.5*dx));


% for u = y_pt:y_pb
%     disp(['u = ' int2str(u)]);
%     for v = 1:(3*n_x)
%         
%         uu = u;
%         vv = v - n_x;
%         
%         R3 = ((vv-X1).^2 + (uu-Y1).^2 ).^(3/2);
%         Ex = Ex + 1/(4*pi*eps) * rho3(u,v) * (X1-vv) ./ R3;
%         
%         R3 = ((vv-X2).^2 + (uu-Y2).^2 ).^(3/2);
%         Ey = Ey + 1/(4*pi*eps) * rho3(u,v) * (Y2-uu) ./ R3;
%         
%         
%     end
% end



% Ex = sigma*(1/dx)^2 * Ex;
% Ey = sigma*(1/dx)^2 * Ey;

% Ex = max(abs(1000*Ey0(:)))/max(abs(Ex(:))) *  Ex;
% Ey = max(abs(1000*Ey0(:)))/max(abs(Ey(:))) *  Ey;


if sig1 > 0 && ~no_slayer
    Ey(1:(y_pt-1),:) = -repmat( fliplr(Ey30(1:(y_pt-1))).', [1 n_x]);
    Ey((y_pb+2):end,:) = repmat( (Ey40(1:size(Ey((y_pb+2):end,:),1))).', [1 n_x]);
    Ey((y_pt+1):(y_pb-0), x_pl:x_pr)...
        = repmat( Ey30(1:length((y_pt+1):(y_pb-0))).', [1 length(x_pl:x_pr)]) ...
        - repmat( flipud(Ey40(1:length((y_pt+1):(y_pb-0))).'), [1 length(x_pl:x_pr)]);
    Ex( y_pt:y_pb, (x_pl+1):x_pr )...
        = grad1 * ( Ey30(1:length((x_pl+1):x_pr))-fliplr(Ey30(1:length((x_pl+1):x_pr))) ) ...
        + grad2 * ( Ey40(1:length((x_pl+1):x_pr))-fliplr(Ey40(1:length((x_pl+1):x_pr))) );
end




%% recation

R = zeros(size(C(:,:,1)));

%% detectors

convg_sim = zeros(1, n_t);
mass_sim = zeros(1, n_t);
chrg_sim = zeros(1, n_t);

c_dtr = zeros(n_i, n_t);

jy_dtr = zeros(n_i, n_t);

%% iterations

% Ex = E_fx;
% Ey = E_fy;

if enable_sym
    x_sym = (size(C,2)+1)/2;
else
    x_sym = size(C,2);
end

y2 = linspace(dx,l_y, size(C,1));
y1 = linspace(0.5*dx,l_y+0.5*dx, size(C,1)+1);
tic;

C0 = C;
C_delta = 0;

filename = [label '_' int0str(cc,4) 'nM_' int0str(1000*c_salt,4) 'mM_k' int0str(kk, 2) '_sig' ...
          int0str(1000*sig1, 4) '_pore' int0str(10*pore_size, 2) ...
          '_y' int2str(amo_y) 'x' int2str(amo_x)];

if real_time_movie>0
    figure(999);
end


if enable_mex
    

    n_mex = ceil( n_t/n_display);
    c_dtr = zeros(n_i, n_mex);
    jy_dtr = zeros(n_i, n_mex);
    c_conv_er = 999;

    for tt = 1:n_mex
        
        a = dt*faraday/eps;
        b = faraday/r/tmp;
        c = kr * dt / k_m / (6e23) / dx^3;
        
        [C, Ex, Ey, Jx, Jy]...
            = mex_pnp( C, Ex, Ey, Jx, Jy, ...
              z, d_m, dx, dt, a, b, c, ...
              [c_bulk_pls; c_salt_min; c_salt_pls], ...
              amo_y, amo_x, [y_pt, y_pb, x_pl, x_pr], ...
              x_sym, n_display, no_slayer);
	c_dtr(:,tt) = squeeze( 1e6 * C(574, 111, :) );
        jy_dtr(:,tt) = squeeze( Jy(400, 111, :) );


	if tt > 1
	    c_conv_er = ( c_dtr(1,tt) - c_dtr(1,tt-1) ) / c_dtr(1,tt-1);
        end

	t = tt*n_display;
	disp(['CONVG = ' num2str(c_conv_er) ' @ step ' int2str(t)]);

        t_now = toc;
        eta = t_now/t * (n_t-t);
        
        t_now = toc;
        eta = t_now/t * (n_t-t);
        disp([num2str(round(1000*t/n_t)/10) '% done, ' num2str(round(eta)) 's to go...']);
        
        if mod(t, 5000000) == 0
            disp('saving intermediate file');
            save([filename '_t' int0str(t, 9)  '.mat'], ...
              'C', 'c_dtr', 'jy_dtr', 'Jx', 'Jy', ...
              '-v7.3');
        end    

        
        
    end
    
    
else
    

    for t = 1:n_t

        C_pv = C;

        % update E 

        if enable_mex
            [Ex, Ey] = mex_exey(Ex, Ey,...
                Jx(:,:,1), Jx(:,:,2), Jx(:,:,3),...
                Jy(:,:,1), Jy(:,:,2), Jy(:,:,3),...
                x_sym, (dt*faraday/eps), z);
        else
            Ex(:,1:x_sym) = Ex(:,1:x_sym) - (dt*faraday/eps)...
                * (  z(1) * Jx(:,1:x_sym,1) + z(2) * Jx(:,1:x_sym,2) + z(3) * Jx(:,1:x_sym,3) );
            Ey(:,1:x_sym) = Ey(:,1:x_sym) - (dt*faraday/eps)...
                * (  z(1) * Jy(:,1:x_sym,1) + z(2) * Jy(:,1:x_sym,2) + z(3) * Jy(:,1:x_sym,3) );
        end


        % update J

        if enable_mex


            [Jx(:,:,1), Jy(:,:,1)] = mex_jxjy(-d_nh4/dx, C(:,:,1), (faraday/r/tmp)*d_nh4*z(1)/2, Ex, Ey, x_sym);
            [Jx(:,:,2), Jy(:,:,2)] = mex_jxjy( -d_cl/dx, C(:,:,2), (faraday/r/tmp)*d_cl*z(2)/2,  Ex, Ey, x_sym);
            [Jx(:,:,3), Jy(:,:,3)] = mex_jxjy( -d_na/dx, C(:,:,3), (faraday/r/tmp)*d_na*z(3)/2,  Ex, Ey, x_sym);

        else

            if enable_sym
                for ii = 1:n_i
                    Jx(:,2:x_sym,ii) = - d_m(ii)/dx * diff(C(:,1:x_sym,ii), 1, 2) ...
                        + 1*(faraday/r/tmp) * d_m(ii) * z(ii) * 0.5 * ( C(:,1:(x_sym-1),ii) + C(:,2:(x_sym),ii) ) .* Ex(:,2:x_sym);  

                    Jy(2:(end-1),1:x_sym,ii) = - d_m(ii)/dx * diff(C(:,1:x_sym,ii), 1, 1)  ...
                        + 1*(faraday/r/tmp) * d_m(ii) * z(ii) * 0.5 * ( C(1:(end-1),1:x_sym,ii) + C(2:end,1:x_sym,ii) ) .* Ey(2:(end-1),1:x_sym);

                end
            else
                Jx(:,:,1) = - d_nh4/dx * diff([C(:,end,1) C(:,:,1) C(:,1,1)], 1, 2) ...
                    + 1*(faraday/r/tmp) * d_nh4 * z(1) * 0.5 * ( [C(:,end,1) C(:,:,1)] + [C(:,:,1) C(:,1,1)] ) .* Ex(:,:);  
                Jx(:,:,2) = - d_cl/dx * diff([C(:,end,2) C(:,:,2) C(:,1,2)], 1, 2) ...
                    + 1*(faraday/r/tmp) * d_cl* z(2) * 0.5 * ( [C(:,end,2) C(:,:,2)] + [C(:,:,2) C(:,1,2)] ) .* Ex(:,:);  
                Jx(:,:,3) = - d_na/dx * diff([C(:,end,3) C(:,:,3) C(:,1,3)], 1, 2) ...
                    + 1*(faraday/r/tmp) * d_na * z(3) * 0.5 * ( [C(:,end,3) C(:,:,3)] + [C(:,:,3) C(:,1,3)] ) .* Ex(:,:);  

                Jy(2:(end-1),:,1) = - d_nh4/dx * diff(C(:,:,1), 1, 1)  ...
                    + 1*(faraday/r/tmp) * d_nh4 * z(1) * 0.5 * ( C(1:(end-1),:,1) + C(2:end,:,1) ) .* Ey(2:(end-1),:);
                Jy(2:(end-1),:,2) = - d_cl/dx * diff(C(:,:,2), 1, 1)  ...
                    + 1*(faraday/r/tmp) * d_cl * z(2) * 0.5 * ( C(1:(end-1),:,2) + C(2:end,:,2) ) .* Ey(2:(end-1),:);
                Jy(2:(end-1),:,3) = - d_na/dx * diff(C(:,:,3), 1, 1)  ...
                    + 1*(faraday/r/tmp) * d_na * z(3) * 0.5 * ( C(1:(end-1),:,3) + C(2:end,:,3) ) .* Ey(2:(end-1),:);
            end

        end

        % Diffusive BCs

        Jy(1,1:x_sym,1) = - d_nh4/dx * ( C(1,1:x_sym,1) - c_bulk_pls(1,1:x_sym,1) - C_delta);
        Jy(1,1:x_sym,2) = - d_cl /dx * ( C(1,1:x_sym,2) - c_salt_min(1,1:x_sym,1));
        Jy(1,1:x_sym,3) = - d_na /dx * ( C(1,1:x_sym,3) - c_salt_pls(1,1:x_sym,1));    


        % Neumann BCs

        Jx(y_pt:y_pb, x_pl+1,:) = 0;
        Jx(y_pt:y_pb, x_pr,:) = 0;

        Jy(y_pt, [1:x_pl x_pr:end],:) = 0;
        Jy(y_pb+1, [1:x_pl x_pr:end],:) = 0;

        Jy(end,:,:) = 0;


        % update C
        if enable_mex

            C(:,:,1) = mex_c( C(:,:,1), -dt/dx, Jx(:,:,1), Jy(:,:,1), x_sym);
            C(:,:,2) = mex_c( C(:,:,2), -dt/dx, Jx(:,:,2), Jy(:,:,2), x_sym);
            C(:,:,3) = mex_c( C(:,:,3), -dt/dx, Jx(:,:,3), Jy(:,:,3), x_sym);

        else
            if enable_sym
                C(:,1:x_sym,:) = C(:,1:x_sym,:) - dt/dx * ( diff([Jx(:,1:x_sym,:) -Jx(:,x_sym,:)], 1, 2)...
                    + diff(Jy(:,1:x_sym,:), 1, 1) );
            else
                C = C - dt/dx * ( diff(Jx, 1, 2) + diff(Jy, 1, 1) );
            end
        end
        R(amo_y,amo_x) =  kr * dt * (C(amo_y,amo_x,1)/k_m) / (6e23) / dx^3;

        C(:,:,1) = C(:,:,1) - R;
        C_delta = 0; % sum(R(:)) / n_x;
        C(C<0) = 0;
    %     C(1,:,:) = c_bulk;

    %     C(1,:,1) = C(1,:,1) - sum(R(:))/n_x;


        % monitors
        c_dtr(:,t) = squeeze( C(amo_y,amo_x,:) );
        if t == 1
            c_conv_er = (c_dtr(1,t) - C0(amo_y, amo_x, 1)) / C0(amo_y, amo_x, 1);
        else
            c_conv_er = diff(c_dtr(1, (t-1):t), 1, 2) / c_dtr(1, t-1);
        end
        convg_sim(t) = c_conv_er;


        jy_dtr(:,t) = squeeze( Jy( round(y_pc/dx), round((n_x+1)/2), :) );

        mass_sim(t) = sum( C(:) );
        chrg_sim(t) = sum( sum( C(:,:,1) - C(:,:,2) ) );


        if mod(t,n_display) == 0
            disp(['CONVG = ' num2str(c_conv_er) ' @ step ' int2str(t)]);
            t_now = toc;
            eta = t_now/t * (n_t-t);

            if c_conv_er == 0
                disp(['early convergence occured at ' int2str(t) 'th step,'...
                    ' or simulation time ' num2str(1e9*t*dt) 'ns.']);
                mass_sim((t+1):end) = mass_sim(t);
                chrg_sim((t+1):end) = chrg_sim(t);
                c_dtr(:,(t+1):end) = c_dtr(:,t);
            break;
            else
                disp([num2str(round(1000*t/n_t)/10) '% done, ' num2str(round(eta)) 's to go...']);
            end

            if real_time_movie>0
                plot(y2, 0.001*C0(:,111,1), 'k', y2, 0.001*C(:,111,1));
    %             ylim([0 2e-9*]);


        %         imagesc( C(:,:,1) );
        %         caxis([0 2]);
            %     colorbar;
                title(['t = ' int2str(t) '/' int2str(n_t)]);
    %             ylim([v_z 0])
        %         colormap(jet);
        %         hold on;
                pause(0.01);
            %     print('-dpng', ['pot_sc050_' int0str(t,4) '.png']);
            end
        end
        if mod(t, 500000) == 0
        disp('saving intermediate file');
            save([filename '_t' int0str(t, 7)  '.mat'], ...
              'C', 'c_dtr', 'jy_dtr', 'Jx', 'Jy', 'Ex', 'Ey',...
              '-v7.3');
        end    



    end
end

% hold off;
toc

if enable_sym
    C = [C(:,1:x_sym,:) fliplr(C(:,1:(x_sym-1),:))];
    Jx = [Jx(:,1:x_sym,:) -fliplr(Jx(:,1:x_sym,:))];
    Jy = [Jy(:,1:x_sym,:) fliplr(Jy(:,1:(x_sym-1),:))];
    Ex = [Ex(:,1:x_sym,:) -fliplr(Ex(:,1:x_sym,:))];
    Ey = [Ey(:,1:x_sym,:) fliplr(Ey(:,1:(x_sym-1),:))];
end



P = -cumsum( dx * Ey(1:end,:) );


save([filename '.mat'],...
    'C0','C','P','Ex','Ey','Jx','Jy',...
    'convg_sim','mass_sim','chrg_sim','c_dtr','jy_dtr','n_display',...
    '-v7.3');
% 

end
