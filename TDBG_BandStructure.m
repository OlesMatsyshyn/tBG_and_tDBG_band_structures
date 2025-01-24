clear all
% This code in many places uses the notation and values from the paper:
% Maximally Localized Wannier Orbitals and the Extended Hubbard Model for
% Twisted Bilayer Graphene
% PHYSICAL REVIEW X 8, 031087
% Open it, and have side by side for the better understanding.

%% Parameters, constants 
theta  = 1.0 * pi / 180; % twist angle in radians
valley = 1;               % valley index
Delta  = 0;               % stagger potential induced by hBN 
alpha  = 2135.4;          % meV, is the hbar*v_F/a; look after Eq.(2) in 
                          % section EFFECTIVE CONTINUUM MODEL

theta_str = join(['\theta = ', num2str(theta*180/pi,'%4.2f')]);

% The displayed band structure is along some path in k space
% in this code its K_+ -> Gamma -> M -> K_-, with 3 intervals
% here set the number of points you want for each interval:
k_size_K_to_Gamma = 100;
k_size_Gamma_to_M = 100;
k_size_M_to_Kprim = 50;

% Select the sutable number. Controlls the accuracy of the model.
% N_moire is a positive integer. Complexity grows fast with
% N_moire. For N_moire = 10, code takes a few minutes to execute.
N_moire = 3; % indicates how many cites (forward and backward) 
             % in k space are coupled at given k_now 
             % resulting H, will have 8*(2*N_moire+1)^2 bands

% Pauli Matrises
sx = [0, 1 ;
      1, 0];
sy = [ 0, -1i ;
      1i,   0];
sz = [1, 0 ;
      0,-1];
s0 = eye(2); 

%% ATOMIC STRUCTURE & the graphene lattice, the Moiré lattice
% We employ the lattice convention from PHYSICAL REVIEW X 8, 031087
% See the section ATOMIC STRUCTURE
%
% as in the paper, lattice vectors in units of graphene lattice constant:
a1 = [1; 0]; a2 = [1/2; sqrt(3)/2]; %  note a = 0.246 nm (irrelevant) 
% reciprocal lattice vectors in units of 1/a [also denoted as b1, b2]
a1_star = 2*pi*[1; -1/sqrt(3)]; a2_star = 2*pi*[0; 2/sqrt(3)]; 

% Wavevector b [see below] points to one of the K points in
% the Brillouin zone of graphene, is aligned along the x-axis
% location of K = -valley*(2*a1 + a2)/3, before the rotation
b_Kx = 4*pi/(3)*[1;0];
% Missmatch of these vectors in top and bottom layer define 
% the moiré reciprocal lattice size:

% Anti-clock wise rotation by angle th in 2D
R2D = @(th)[ cos(th), -sin(th)  ;
             sin(th),  cos(th) ];
% Create the rotational matrises with values
R1 = R2D(-theta/2); % Rotation of vectors in layer 1
R2 = R2D(+theta/2); % Rotation of vectors in layer 2
% Thus we find the missmath as
b_M = norm(R1*b_Kx - R2*b_Kx);

% Rotated packed vectors are
% access by ( layer (1 or 2) , component (x or y), vector number (1 or 2))  
a_th = zeros(      2         ,       2           ,            2          );
a_th(1,:,:) = [R1*a1,R1*a2];
a_th(2,:,:) = [R2*a1,R2*a2];
% thus a_th(1,:,1) = rotated a1 in layer 1;
%      a_th(2,:,2) = rotated a2 in layer 2
% Similarly
a_th_star = zeros(2,2,2); 
a_th_star(1,:,:) = [R1*a1_star,R1*a2_star];
a_th_star(2,:,:) = [R2*a1_star,R2*a2_star];

% Moiré vectors, using the definitions directly from the text
G1M_paper = a_th_star(1,:,1) - a_th_star(2,:,1);
G2M_paper = a_th_star(1,:,2) - a_th_star(2,:,2);
% Or just use the known answer (as we do :) ):
G1M = [-1/2;-sqrt(3)/2]*sqrt(3)*b_M;
G2M = [   1;         0]*sqrt(3)*b_M;
% these are actual values for the missmatched K points
% but they are irrelevant for the moire BZ.
K1_G1 = -valley*R1*b_Kx; % Graphene top layer K point
K2_G1 = -valley*R2*b_Kx; % Graphene bot layer K point
% Shifted by q0 = (K1_G1 + K2_G1)/2; K points, see text after Eq.(A4)
K1_paper = K1_G1 - (K1_G1 + K2_G1)/2 + G2M_paper'/2;
K2_paper = K2_G1 - (K1_G1 + K2_G1)/2 + G2M_paper'/2;

K1_moire  = valley*[sqrt(3)/2; 1/2]*b_M;
K2_moire  = valley*[sqrt(3)/2;-1/2]*b_M;
%% The path in the k-space is      : K_+ -> Gamma ->    M      -> K_-
%  PHYSICAL REVIEW X 8, 031087
%                         notation : K   -> Gamma ->    M      -> K^\prime
%          or in our code notation : K1  -> [0 0] -> (K1+K2)/2 -> K2        
%  See Fig 1b insert of PHYSICAL REVIEW X 8, 031087

% Points on the effective (dimentionless) moiré Brilluin zone
%          _K2_
%     _M''/    \M' _
%   _/              \_
% K1                  K1
% |                    |
% M          G         M 
% |                    |
% K2_                _K2
%    \_ M'      M''_/
%         \_  _/
%           K1
K1_eff    = [sqrt(3)/2;1/2];
K2_eff    = [sqrt(3)/2;-1/2];
Gamma_eff = [0 0];
M_eff     = (K1_eff + K2_eff)/2;

% number of k points in each peace of the path. Tunable
% number of k points in each peace of the path. Tunable
N_path = [k_size_K_to_Gamma, k_size_Gamma_to_M, k_size_M_to_Kprim];
% the path is empty, at first
k_grid = []; 
% now we add whatever trajectory we want
k_grid = path_add(k_grid,K1_eff   ,Gamma_eff,N_path(1));
k_grid = path_add(k_grid,Gamma_eff,M_eff    ,N_path(2));
k_grid = path_add(k_grid,M_eff    ,K2_eff   ,N_path(3));
k_length = length(k_grid);

%% EFFECTIVE CONTINUUM MODEL
% We employ the convention and values from PHYSICAL REVIEW X 8, 031087
% See section EFFECTIVE CONTINUUM MODEL
%
% note u1 = u, u2 = u^\prime
% with interlayer tunnellings:
u1  = 79.7;   % meV
u2  = 97.5;   % meV
% gammas
gamma1 = 400; % meV 400
gamma3 = 320; % meV 320
gamma4 = 44;  % meV
d_prime = 50; % meV
d = 20;       % meV

v3 = (sqrt(3)/2)*gamma3;
v4 = (sqrt(3)/2)*gamma4; 

S = @(k) [v4 * (k(1)*valley + 1j * k(2)), gamma1;
          v3 * (k(1)*valley - 1j * k(2)), v4 * (k(1)*valley + 1j * k(2))];
H0       = @(k) -alpha * (k(1)*valley * sx + k(2) * sy) + ...
                                                 [ 0 , 0       ;
                                                   0 , d_prime];
H0_prime = @(k) -alpha * (k(1)*valley * sx + k(2) * sy) + ...
                                                 [ d_prime , 0 ;
                                                   0       , 0];
% We inititalise the tunnelling matrix as in Eq(3) of the paper.
% The position dependent factors are irrelevet here, as they dictate 
% how we fill up the full Hamiltonian, with these constant matrises:
fai = valley*2*pi/3;
T1 = [u1, u2 ;
      u2, u1];
 
T2 = [u1,             u2*exp(-1i*fai);
      u2*exp(1i*fai), u1            ];

T3 = [u1,              u2*exp(1i*fai);
      u2*exp(-1i*fai), u1           ];

% Constructing the Hamiltonian
% Number of moiré lattice jumps back and forward to consider

m = generate_moire_lattice_indexing(N_moire);
N_insert = (2*N_moire+1)^2; % Number of matrix insertions
                            % 2*Nbands is the number of bands in one layer
N_bands = 8*N_insert;

Energy = zeros(N_bands, k_length);

% Further matrix sizes are exclusively defined by the value N_moire
% essentially  for every k as a center [0], we will build an
% effective Hamiltonian coupled to the k's shifted by
% moire lattice vectors G1, G2. More shifts we include - better the 
% model. For more info Read section EFFECTIVE CONTINUUM MODEL. 
% Example: for N_moire = 1 we couple momentums in the following way:
% -G2 +0 +G2
%  o-->o-->o    <- +G1
%    \   \
%  o-->0-->o    <- + 0
%    \   \  
%  o-->o-->o    <- -G1
%
% You can note, that to include all the horisontal couplings, 
% arrows just point to one direction (e.g. left). 
% Do not double include those! 
% Hamiltonians at different k's couple either 
%       1. On cite                      (T1). 
%       2. Shifted by + G1              (T2).
%       3. Shifted by (both!) + G1 + G2 (T3).
% [see Eq.(3) of EFFECTIVE CONTINUUM MODEL]: 

textprogressbar('>> calculating the band structure: ');
% For each point in K space we want, we need to solve the effective H
for iK = 1 : k_length
   textprogressbar(fix(100*iK/k_length));
   % Current value of k. k-scale is defined by b_M
   k_now = k_grid(:,iK)*b_M; 
   % Matrix place holders
   U  = zeros(4*N_insert,4*N_insert);
   H1 = zeros(4*N_insert,4*N_insert);
   H2 = zeros(4*N_insert,4*N_insert);
   % The idea is to build the following Hamiltonian for a given k 
   % where U(k1,k2), is how k1 is coupled to k2
   % H1[kt1] S'[kt1]   0      0     ..0        0      0       0.......
   % S[kt1]  H2[kt1]   0      0     ..0        T0     0       0....... 
   % 0       0       H1[kt2] S'[kt2]..0        0      0       0.......
   % 0       0       S[kt2]  H2[kt2]..0        T1     0       0.......
   % .................................0        0      0       0.......
   % .................................................................
   % --------------------------------H3[kb1] S'[kb1]   0       0...... 
   % ---------U^\dagger--------------S[kb1]  H4[kb1]   0       0......
   % --------------------------------0           0   H3[kb2] S'[kb2]..
   % --------------------------------0           0   S[kb2]  H4[kb2]..
   % .................................................................
   % with U[k1,k1] = T1; U[k1,k1 + delta_k] = T2; U[k1,k1 + 2*delta_k] = T3

   for i = 1 : N_insert

       % We define the momentum as in Eq.(2) and see the paragraph before 
       % Eq.(4) of the section "EFFECTIVE CONTINUUM MODEL" 
       % k in layer top, rotation of k vectors is opposite to a
       kl1 = R2 * (k_now - K1_moire + m(i,1) * G1M + m(i,2) * G2M); 
       % k in layer bot
       kl2 = R1 * (k_now - K2_moire + m(i,1) * G1M + m(i,2) * G2M); 
       
       % Diagonal components of the blocks
       H1(4*i-3:4*i,4*i-3:4*i) = [ H0(kl1)+3*s0*d/2, S(kl1)'              ;
                                   S(kl1)          , H0_prime(kl1)+s0*d/2];
       H2(4*i-3:4*i,4*i-3:4*i) = [ H0(kl2)-s0*d/2, S(kl2)'                ;
                                   S(kl2)        , H0_prime(kl2)-3*s0*d/2]; 
       U(4*i-3:4*i,4*i-3:4*i)  = [ zeros(2,2), T1         ;
                                   zeros(2,2), zeros(2,2)];  
       for j = 1 : N_insert
           % horisontal arrows, for xi = 1, x_i = x_j + 1, y_i = y_j
           if m(j,1)-m(i,1)==-valley*1 && m(i,2)==m(j,2)
               U(4*i-3:4*i,4*j-3:4*j) = [ zeros(2,2), T2         ;
                                          zeros(2,2), zeros(2,2)];
           end
           if m(j,1)-m(i,1)==-valley*1 && m(j,2)-m(i,2)==-valley*1
               U(4*i-3:4*i,4*j-3:4*j) = [ zeros(2,2), T3         ;
                                          zeros(2,2), zeros(2,2)];
           end     
       end
   end    
   Energy(:,iK) = eig([H1 U' ;
                       U H2]);
end
textprogressbar(' done. plotting...');

% displaying the results
figure
set(gcf,'Position',[100 300 300 500])
for i = 1:height(Energy)
    hold on;
    plot(Energy(i,:),'-','Linewidth',1.1,'Color','b');
    ylabel('E(meV)');
    set(gca,'xticklabel',[])
    box on
end
title(theta_str)
axis([0 301 -80 100]);
set(gca,'XTick',[0, k_size_K_to_Gamma,...
                  k_size_Gamma_to_M+k_size_K_to_Gamma, k_length]);
xlim([0 k_length])
ylim([-50 100])
set(gca,'YTick',-300:50:300);
set(gca,'XTickLabel',...
   {'$\rm\bar{K}$';'$\bar{\Gamma}$';'$\rm\bar{M}$';'$\rm\bar{K}^\prime$'}); 
set(gca,'TickLabelInterpreter','latex')

hold on 
plot([100 100],[-300 300],'--','LineWidth',0.9,'Color','k');
hold on
plot([200 200],[-300 300],'--','LineWidth',0.9,'Color','k');

% Generate moiré reciprocal lattice points indexing

function lattice_points = generate_moire_lattice_indexing(N)
    lattice_points = zeros((2*N + 1)^2, 2);
    index = 1;
    for i = -N:N
        for j = -N:N
            lattice_points(index, :) = [i, j];
            index = index + 1;
        end
    end
end

function path = path_add(A,p1,p2,N)
    if isempty(A)
        path = [linspace(p1(1),p2(1),N+1); linspace(p1(2),p2(2),N+1)];
        return 
    end
    x = linspace(p1(1),p2(1),N+1);
    y = linspace(p1(2),p2(2),N+1);
    path = [A(1,:) x(2:end); A(2,:) y(2:end)];
end

% Just a fancy progressbar. Delete/cpmment it if complicated.

function textprogressbar(c)
    % Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
    % Version: 1.0
    % Changes tracker:  29.06.2010  - First version
    %% Initialization
    persistent strCR;          % Carriage return pesistent variable
    
    % Vizualization parameters
    strPercentageLength = 10;  % Length of percentage string (must be >5)
    strDotsMaximum      = 10;  % The total number of dots in a progress bar
    
    %% Main 
    
    if isempty(strCR) && ~ischar(c)
        % Progress bar must be initialized with a string
        error('The text progress must be initialized with a string');
    elseif isempty(strCR) && ischar(c)
        % Progress bar - initialization
        fprintf('%s',c);
        strCR = -1;
    elseif ~isempty(strCR) && ischar(c)
        % Progress bar  - termination
        strCR = [];  
        fprintf([c '\n']);
    elseif isnumeric(c)
        % Progress bar - normal progress
        c = floor(c);
        percentageOut = [num2str(c) '%%'];
        percentageOut = [percentageOut repmat(' ',1,...
                             strPercentageLength-length(percentageOut)-1)];
        nDots = floor(c/100*strDotsMaximum);
        dotOut = ['[' repmat('.',1,nDots) ...
                                   repmat(' ',1,strDotsMaximum-nDots) ']'];
        strOut = [percentageOut dotOut];
        
        % Print it on the screen
        if strCR == -1
            % Don't do carriage return during first run
            fprintf(strOut);
        else
            % Do it during all the other runs
            fprintf([strCR strOut]);
        end
        
        % Update carriage return
        strCR = repmat('\b',1,length(strOut)-1);
        
    else
        % Any other unexpected input
        error('Unsupported argument type');
    end
end

