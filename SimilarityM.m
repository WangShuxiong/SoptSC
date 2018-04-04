function W = SimilarityM(X,lambda,data)
% Solving the following optimization problem by ADMM
%
%    min_{Z,E}  ||Z||_* + lambda ||E||_{2,1}
%    s.t.       X = XZ + E;
%               Z'1 = 1;
%               Z_{i,j} = 0 for (i,j)\in Omega
%
% Input
%   1) Single cell data X, a m*n matrix, with m rows(genes) and n columns(cells);
%   2) lambda, the default value is 0.5;
%   3) K: number of nearest neighbor points, which is set as
%      K = max(10,min(0.2*m,20))
%
% Output
%   1) Cell-to-cell similarity matrix for clustering: W
%   2) Cell-to-cell similarity matrix for pseudotime: P


[m,n] = size(X);
if m>=60
    [coeff1,X1,pca_eigvalue1] = pca(X','NumComponents',60);
else
    [coeff1,X1,pca_eigvalue1] = pca(X','NumComponents',m);
end

% covx = cov(X');
% pca_eigvalue1 = sort(eig(covx),'descend');
[~,No_Comps1] = max(abs(pca_eigvalue1(2:end-1) - pca_eigvalue1(3:end)));

display(sum(pca_eigvalue1(1:No_Comps1+1))./sum(pca_eigvalue1))

display(No_Comps1);

cc = cumsum(pca_eigvalue1(2:end));
dd = cc(2:end)./sum(pca_eigvalue1(2:end));

K1 = length(find(dd<=0.3));
% display(K1);

if K1 <= 10
    K = 10;
elseif K1 >=30
    K = 30;
else
    K = K1+1;
end
InitY = pca1(data',3);
X2 = tsne(X','Standardize',true,'Perplexity',20,'NumDimensions',3,'InitialY',InitY);
% X2 = tsne(X', [], InitY, [], 20);

display(K);

D = ones(n,n);
if No_Comps1>=1
    No_Comps1 = 1;
end

No_Comps1 = 1;
[IDX,~] = knnsearch(X2(:,1:No_Comps1+2),X2(:,1:No_Comps1+2),'k',K);

for jj = 1:n
    D(jj,IDX(jj,:)) = 0;
end

Z = computM(D,X,lambda);
Z(Z<=eps) = 0;
W = 0.5.*(abs(Z)+abs(Z'));


    function Z = computM(D,X,lambda)
        %% ADMM iteration
        maxiter = 100;
        Err = zeros(maxiter,2);
        rho = 5;        % 5
        mu = 10^(-6);   % 10^(-6)
        mumax = 10^(6);
        epsilon = 10^(-5);
        [m,n] = size(X);
        Z = zeros(n); E = zeros(m,n);
        Y1 = zeros(m,n); Y2 = zeros(1,n); Y3 = zeros(n,n);
        iter = 0;
        
        while 1
            iter = iter + 1;
            if iter >= maxiter
                break;
            end
            
            % step 1: Update J
            mu = min(rho*mu,mumax);
            [U,S,V] = svd(Z-Y3/mu);
            
            R = length(diag(S));
            Dmu = zeros(R,R);
            MM = max(diag(S)-(1/mu)*ones(R,1),zeros(R,1));
            for i = 1:R
                Dmu(i,i) = MM(i);
            end
            
            J = U*Dmu*V';
            
            if iter >= 3
                if Err(iter-1,1) >= Err(iter-2,1) || norm(X-X*Z) <= epsilon
                    break;
                end
            end
            
            
            % Update E
            Q = X-X*Z+Y1/mu;
            for j = 1:n
                if norm(Q(:,j)) > lambda/mu
                    E(:,j) = ((norm(Q(:,j))-lambda/mu)/norm(Q(:,j)))*Q(:,j);
                else
                    E(:,j) = zeros(length(Q(:,j)),1);
                end
            end
            
            eta = norm(X)^2 + norm(ones(n,1))^2+1;
            H = -X'*(X - X*Z - E + (1/mu)*Y1) - ones(n,1)*(ones(n,1)'- ones(n,1)'*Z + (1/mu)*Y2) + ...
                (Z - J + (1/mu)*Y3);
            Z = Z - (1/eta)*H;
            Z(D>0) = 0;
            
            % Update Dual variable
            Y1 = Y1 + mu*(X-X*Z-E);
            Y2 = Y2 + mu*(ones(1,n) - ones(1,n)*Z);
            Y3 = Y3 + mu*(Z - J);
            
            fprintf('%d, %8.6f, %8.6f\n',iter,norm(X-X*Z-E),norm(Z-J));
            Err(iter,:) = [norm(X-X*Z-E) norm(Z-J)];
            
            if max(Err(iter,:)) <= epsilon
                break;
            end
        end
    end
end