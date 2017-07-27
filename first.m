clear;
M = 16; %no. of users
G = 3; %no. of user groups
K = G;
%g{1} = [1,4,6,9,16]; g{2} = [2,7,8,11,14]; g{3} = [3,5,10,12,13,15];
g{1} = [1,2,3,4,5]; g{2} = [6,7,8,9,10]; g{3} = [11,12,13,14,15,16];
N = 6; %no. of transmit antennas
theta = (360/M):(360/M):360; %angular distribution of users (mean)
ang_spread = 1:15; %spread angle of each user
sigma = 1; %noise power
gamma = .01; %QoS

for a = 1:length(ang_spread)
    
    for e = 1:2
        R = zeros(N,N,M); %channel covariance matrix
        eig_val = cell(M,1);
        eig_vec = cell(M,1);
        N_set_bar = cell(M,1);
        N_bar = zeros(M,1);
        for m = 1:M
            for c = 1:N
                for r = 1:N
                    R(r,c,m) = exp(1i*pi*(r-c)*sin(theta(m)))*...
                        exp(-((pi*(r-c)*ang_spread(a)*cos(theta(m)))^2 / 2));
                end
            end
            [eig_vec{m},eig_val{m}] = eig(R(:,:,m));
            eig_val{m} = diag(eig_val{m});
            %N_set_bar{m} = 1:length(eig_val{m});
            N_bar(m) = length(eig_val{m});
        end
        H = cell(6,3);
        Ii = eye(G);
        for k = 1:K
            for m = 1:M
                if (ismember(m,g{k}))
                    A = [zeros(1,N*G); kron(Ii-diag(Ii(:,k)),R(:,:,m))];
                    H{m,k} = [zeros((N*G)+1,1), A];
                    H{m,k}(1,1) = sigma;
                end
            end
        end
        [obj_val,z,eig_val,eig_vec] = feasSearch(g, M, H, eig_val, eig_vec, N_bar, gamma);
        if (obj_val > 0)
            [Mxm, MxI] = max(z);
            for k = 1:K
                if (ismember(MxI,g{k}))
                    ri = find(g{k} == MxI);
                    g{k}(ri) = [];
                    M = M - 1;
                    break
                end
            end
        elseif (obj_val == 0)
            [w_new,t_new] = powerOptimal(H, eig_val, eig_vec, N_bar, gamma);
        end
        continue
    end
        
    
    pow(a) = 20*log10(t_new);
end
plot(ang_spread,pow,'-*r');
