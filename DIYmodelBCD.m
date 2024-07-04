% System initialization
%% Define the strained hoppings and the parameters of the Hamiltonian
t1 = 1;       % In the eV units
t2 = 0.3;       % In the eV units
alpha_g = 0.3;% In the eV units
mu = -1;      % chemical potential
U = 0.75;     % Cooper pairs coupling strength. Yanase has 1.5 of gamma_0 for graphene.                       

lE = 0.3;

% U      |      Tc
% --------------------
% 1.50   |     0.400
% 1.20   |     0.250
% 0.90   |     0.110
% 0.75   |     0.035

hx = 0;
hy = 0;
hz = 0;
Temp = 0.1;% The system temperature

fFD  = @(x) 1/2-tanh(0.5*(x)/Temp)/2; % Fermi-Dirac function
dfFD = @(x) -(0.25).*sech(0.5*(x)/Temp).^2; % We should later divide by Temp

% Grid density
nx = 409;
ny = 400;
q = 0.6;
% The f of the strained graphene
xik  = @(kx,ky) -2*t1*(cos(kx)+cos(3*ky))+4*t2*cos(kx).*cos(ky)-mu;
gkx  = @(kx,ky) -alpha_g*sin(ky);
gky  = @(kx,ky)  alpha_g*sin(kx);
gkz  = @(kx,ky) lE * cos(kx).*cos(ky);

% Derivative w.r.t. x
dxxik  = @(kx,ky) 2*t1*sin(kx)-4*t2*sin(kx).*cos(ky);
dxgkx  = @(kx,ky) 0;
dxgky  = @(kx,ky) alpha_g*cos(kx);
dxgkz  = @(kx,ky) -lE * sin(kx).*cos(ky);
% Derivative w.r.t. y
dyxik  = @(kx,ky) 3*2*t1*sin(3*ky)-4*t2*cos(kx).*sin(ky);
dygkx  = @(kx,ky) -alpha_g*cos(ky);
dygky  = @(kx,ky) 0;
dygkz  = @(kx,ky) -lE * cos(kx).*sin(ky);
% Pauli matrises
s0 = [1,  0 ;
      0 , 1];
sx = [0,  1 ;
      1,  0];
sy = [0 ,-1j;
      1j, 0];
sz = [1,  0;
      0 ,-1];

H_N   = @(kx,ky) xik(kx,ky)*s0 + (gkx(kx,ky)-hx)*sx + (gky(kx,ky)-hy)*sy + (gkz(kx,ky)-hz)*sz;
H_NT  = @(kx,ky) xik(kx,ky)*s0 + (gkx(kx,ky)-hx)*sx - (gky(kx,ky)-hy)*sy + (gkz(kx,ky)-hz)*sz;

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


a1 = [2*pi, 0];
a2 = [0, 2*pi];

v_table = zeros (nebands, nebands, 2); 
energy  = zeros (nebands, nx+1, ny+1);
BCD     = zeros(2,2,2);
BCDD    = zeros(nx+1,ny+1,2,2,2);
for i = 0 : nx
    for j = 0 : ny
        k = [-pi,-pi] + (i/nx) * a1 + (j/ny) * a2;
        [Vec,en] = eigenshuffle(H_N(k(1),k(2)));
        energy(:,i+1,j+1) = en(:); % band energy at k
        for n = 1 : nebands
            for m = 1 : nebands
                v_table(n,m,1) = Vec(:,n)'*dxH_N(k(1),k(2))*Vec(:,m);
                v_table(n,m,2) = Vec(:,n)'*dyH_N(k(1),k(2))*Vec(:,m);
                if n == m
                    continue
                end
                for alpha = 1:2
                    for beta = 1:2
                        for gamma = 1:2
                            BCD(alpha,beta,gamma) = BCD(alpha,beta,gamma)+...
                                + (v_table(n,n,alpha)*dfFD(en(n))-v_table(m,m,alpha)*dfFD(en(m)))*...
                                v_table(n,m,beta)*v_table(m,n,gamma)/(en(n)-en(m))^2;

                            BCDD(i+1,j+1,alpha,beta,gamma) = BCDD(i+1,j+1,alpha,beta,gamma)+...
                                + (v_table(n,n,alpha)*dfFD(en(n))-v_table(m,m,alpha)*dfFD(en(m)))*...
                                v_table(n,m,beta)*v_table(m,n,gamma)/(en(n)-en(m))^2;
                        end
                    end
                end
            end
        end
    end
end
dkxdky = 1/(nx*ny); % dkx dky / (2*pi)^2
BCD = BCD * dkxdky/Temp;

surf(squeeze(energy(1,:,:)))
hold on
surf(squeeze(energy(1,:,:)))
% Dq = qx*0+0.1; % seed value
% Jq = qx*0;
% nqx = 30;
% qx = linspace(-0.1,0.1,nqx);
%%  Not yet required 
% tic
% cprintf('hyper','Starting the canlculation... Hold on... It"s not so simple...\n')% XX component (just a nice dislay of the statement)
% parfor iqx = 1:nqx
%     iter = 0;             % number of iterations from this q completed
%     discr = 1;            % just to inititalize while
%     % For the current q find the Bogolubov Band structure:
%     band   = zeros (nbands, nbands, nx+1, ny+1); 
%     energy = zeros (nbands, nx+1, ny+1);
%     while discr > 0.01 && ~isnan(discr) && abs(Dq(iqx)) > 10^(-6) % convergence cryterium
%         iter = iter + 1;  % update iteration index
%         % Do not count the boundary twise! start from 0 and end at nx-1,
%         % or from 1 and end at nx
%         for i = 0 : nx
%             for j = 0 : ny
%                 k = [-pi,-pi] + (i/nx) * a1 + (j/ny) * a2;
%                 [Vec,Val] = eigenshuffle(H_BdG(k(1),k(2),qx(iqx),Dq(iqx)));
%                 band(:,:,i+1,j+1) = Vec;       % column eigenvectors at which kx, at which ky, ;
%                                            % eigenvector vi = band(:,i,kx,ky);
%                 energy(:,i+1,j+1) = Val(:);    % band, at which kx, at which ky;
%             end
%         end
% 
%         Dnm  = @(n,m,ikx,jky) band(:,n,ikx,jky)'*dDH_BdG*band(:,m,ikx,jky);
%                
%         Dmem = Dq(iqx);                 % memorise the D before the update
% 
%         IntegralN = zeros(1,nbands);    % split the task between cores
%         for n = 1 : nbands
%             for i = 0 : nx
%                 for j  = 0 : ny   
%                     IntegralN(n) = IntegralN(n) + fFD(energy(n,i+1,j+1))*Dnm(n,n,i+1,j+1);
%                 end
%             end
%         end
%         % Updating the gap
%         Dq(iqx) = - 0.5 * U * dkxdky * sum(IntegralN);
%         discr = abs((Dq(iqx)-Dmem)/Dmem);
%         %disp(['iteration = ', num2str(iter), ' finished, Delta = ', num2str(fix(Dq(iqx)*1000)/1000),...
%         %                                        ', discrepancy = ',num2str(fix(discr*10^6)/10^6)]);      
%     end
%     if abs(Dq(iqx)) < 10^(-6)
%         Dq(iqx) = 0;
%     end
% %     for i = 0 : nx
% %         for j = 0 : ny
% %             k = [-pi,-pi] + (i/nx) * a1 + (j/ny) * a2;
% %             [Vec,Val] = eigenshuffle(H_BdG(k(1),k(2),qx(iqx),Dq(iqx)));
% %             band(:,:,i+1,j+1) = Vec;       % column eigenvectors at which kx, at which ky, ;
% %                                        % eigenvector vi = band(:,i,kx,ky);
% %             energy(:,i+1,j+1) = Val(:);    % band, at which kx, at which ky;
% %         end
% %     end
% %     % For a converged iterative process compute the current
% %     % Note, the Hamiltonian is k dependent, and the vectors are stored by
% %     % indeces, that is why we need so many entries
% %     Jnm  = @(n,m,ikx,jky,kx,ky,q) band(:,n,ikx,jky)'*dxH_BdG(kx,ky,q)*band(:,m,ikx,jky);
% %     for n = 1 : nbands
% %         for i = 0 : nx
% %             for j  = 0 : ny
% %                 k = [-pi,-pi] + (i/nx) * a1 + (j/ny) * a2;
% %                 Jq(iqx) = Jq(iqx) + fFD(energy(n,i+1,j+1))*Jnm(n,n,i+1,j+1,k(1),k(2),qx(iqx));
% %             end
% %         end
% %     end
%     % displaying progress
%     % progress = fix(100 * iqx / nqx);
%     % dispstat([sprintf('%d%% completed, time passed: ',progress) num2str(fix(10*t)/10)...
%     %            ' sec, time left (approx): ' num2str(fix(10*t*(nqx-iqx)/iqx)/10) ' sec'],progress);
%     disp(['I finished ', num2str(iqx)])
% end
% Jq = Jq * dkxdky;
% plot(qx,Dq)
% xlabel('q_x')
% ylabel('\Delta_q')
% 
% figure
% plot(qx,real(Jq))
% xlabel('q_x')
% ylabel('J_q')
