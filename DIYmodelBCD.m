% System initialization
%% Define the strained hoppings and the parameters of the Hamiltonian
t1x = 1;       % In the eV units
t1y = 1;       % In the eV units
alpha_g = 0.3;% In the eV units
beta_D  = -0.2;
delta   = 0.25;
mu      = -1;      % chemical potential
U       = 0.75;     % Cooper pairs coupling strength. Yanase has 1.5 of gamma_0 for graphene.                       
lE      = 0.1; % to split the bands

% U      |      Tc
% --------------------
% 1.50   |     0.400
% 1.20   |     0.250
% 0.90   |     0.110
% 0.75   |     0.035

Temp = 0.01;% The system temperature
fFD  = @(x) 1/2-tanh(0.5*(x)/Temp)/2; % Fermi-Dirac function
dfFD = @(x) -(0.25).*sech(0.5*(x)/Temp).^2; % We should later divide by Temp

% Grid density, lattice and reciprocal vectors
nx = 409;
ny = 400;
a1 = [2*pi, 0];
a2 = [0, 2*pi];

% The f of the strained graphene
xik  = @(kx,ky) -2*(t1x*cos(kx)+t1y*cos(ky))-mu;
gkx  = @(kx,ky) -alpha_g*sin(ky) - beta_D*sin(kx);
gky  = @(kx,ky)  alpha_g*sin(kx) + beta_D*sin(ky);
gkz  = @(kx,ky) lE * cos(kx).*cos(ky)+delta*(sin(kx)-sin(ky));

% Derivative w.r.t. x
dxxik  = @(kx,ky) 2*t1x*sin(kx);
dxgkx  = @(kx,ky) -beta_D*cos(kx);
dxgky  = @(kx,ky)  alpha_g*cos(kx);
dxgkz  = @(kx,ky) -lE * sin(kx).*cos(ky)+delta*cos(kx);
% Derivative w.r.t. y
dyxik  = @(kx,ky) 2*t1y*sin(ky);
dygkx  = @(kx,ky) -alpha_g*cos(ky);
dygky  = @(kx,ky) beta_D*cos(ky);
dygkz  = @(kx,ky) -lE * sin(ky).*cos(kx)-delta*cos(ky);
% Pauli matrises
s0 = [1,  0 ;
      0 , 1];
sx = [0,  1 ;
      1,  0];
sy = [0 ,-1j;
      1j, 0];
sz = [1,  0;
      0 ,-1];

H_N   = @(kx,ky) xik(kx,ky)*s0 + gkx(kx,ky)*sx + gky(kx,ky)*sy + gkz(kx,ky)*sz;
H_NT  = @(kx,ky) xik(kx,ky)*s0 + gkx(kx,ky)*sx - gky(kx,ky)*sy + gkz(kx,ky)*sz;

dxH_N = @(kx,ky) dxxik(kx,ky)*s0 + dxgkx(kx,ky)*sx + dxgky(kx,ky)*sy + dxgkz(kx,ky)*sz;
dyH_N = @(kx,ky) dyxik(kx,ky)*s0 + dygkx(kx,ky)*sx + dygky(kx,ky)*sy + dygkz(kx,ky)*sz;

H_BdG = @(kx,ky,q,Delta)[H_N(kx+q,ky)      , Delta*1j*sy    ;
                         -conj(Delta)*1j*sy, -H_NT(-kx,-ky)];

dxH_BdG = @(kx,ky,q)[dxH_N(kx+q,ky), zeros(2,2)  ;
                         zeros(2,2), zeros(2,2) ];

dDH_BdG= [zeros(2,2), 1j*sy ;
          -1j*sy    , zeros(2,2)];

nbands = length(dDH_BdG);
nebands = fix(nbands/2);





v_table = zeros (nebands, nebands, 2); 
energy  = zeros (nx+1, ny+1, nebands);
BCD     = zeros(2,2,2);
for i = 0 : nx
    for j = 0 : ny
        k = [-pi,-pi] + (i/nx) * a1 + (j/ny) * a2;
        [Vec,en] = eigenshuffle(H_N(k(1),k(2)));
        % test that energy(i,j+1,:) = energy(nx+1,j+1,:)?

        energy(i+1,j+1,:) = en(:); % band energy at k
        for n = 1 : nebands
            for m = 1 : nebands
                v_table(n,m,1) = Vec(:,n)'*dxH_N(k(1),k(2))*Vec(:,m);
                v_table(n,m,2) = Vec(:,n)'*dyH_N(k(1),k(2))*Vec(:,m);
            end
        end
        for n = 1 : nebands
            for m = 1 : nebands
                if n == m
                    continue
                end
                for alpha = 1:2
                    for beta = 1:2
                        for gamma = 1:2
                            BCD(alpha,beta,gamma) = BCD(alpha,beta,gamma)+...
                                + v_table(n,n,alpha)*dfFD(en(n))*...
                                (v_table(n,m,beta)*v_table(m,n,gamma)-v_table(m,n,beta)*v_table(n,m,gamma))/(en(n)-en(m))^2;

                        end
                    end
                end
            end
        end
    end
end
dkxdky = 1/(nx*ny); % dkx dky / (2*pi)^2
BCD = BCD * dkxdky/Temp;
figure
surf(energy(:,:,1),'edgecolor','none')
hold on
surf(energy(:,:,2),'edgecolor','none')
figure
hold on
plot(energy(fix(end/2),:,1))
plot(energy(fix(end/2),:,2))
disp(['BCD estimation is: ', num2str(sum(abs(BCD),'all')),'; if more than 10e-10, its finite'])

tic
iter = 0;             % number of iterations from this q completed
%% For the current q find the Bogolubov Band structure:
band   = zeros (nbands, nbands, nx+1, ny+1); 
energy = zeros (nbands, nx+1, ny+1);
Dnn    = zeros (nbands, nx+1, ny+1);
Dq = 0.1; % seed value
discr = 1;
while discr > 0.01 && ~isnan(discr) && abs(Dq) > 10^(-6) % convergence cryterium
    iter = iter + 1;  % update iteration index
    % or from 1 and end at nx
    while discr > 10^(-6) && abs(Dq) > 10^(-5) % convergence cryterium
        for i = 0 : nx
            for j = 0 : ny
                k = [-pi,-pi] + (i/nx) * a1 + (j/ny) * a2;
                [Vec,Val] = eigenshuffle(H_BdG(k(1),k(2),0,Dq));
                energy(:,i+1,j+1) = Val(:);    % band, at which kx, at which ky;
                for n = 1 : nbands
                    Dnn(n,i+1,j+1)  = Vec(:,n)'*dDH_BdG*Vec(:,n);
                end
            end
        end     
        Dmem = Dq;                 % memorise the D before the update
        % Updating the gap
        Dq = - 0.5 * U * dkxdky * sum(sum(sum(fFD(energy).*Dnn)));
        discr = abs((Dq-Dmem));
        %disp(['iteration = ', num2str(iter), ' finished, Delta = ', num2str(fix(Dq(iqx)*1000)/1000),...
        %                                        ', discrepancy = ',num2str(fix(discr*10^6)/10^6)]);      
    end     
end
if abs(Dq) < 10^(-6)
    Dq = 0;
end
disp(['Del estimation is: ', num2str(Dq),'; if more than 10e-4, its finite, also should be less then 1'])
