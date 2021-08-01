function [Lambda,Q] = stepwiseCFCPC(X,n,pmax,lmax)
%CPCSTEPWISE calculate the common pricinple components
%Input  X= a cell variable with two items in cell containing Source x
%       frequencies as one dimension and subjects as second dimension      
%       n number of samples for each type of modality or group
%       pmax number of common principle components
%       lmax maximum number of iterations for convergance <p
%Output
%       Lambda   eigenvalues pXk
%       Q              eigenvectors pXp (CPC)
% Reference: Stepwise Common Principal Components Trendafilov

% Usama,Fuleah,Pedro, Andy
[ns,p,k]=size(X);
nt = sum(n);
Q = [];
Wbar = [];
for i=1:k
    Wi(:,:,i) = (X(:,:,i)-mean(X(:,:,i)))/sqrt(n(i));
    Wbar  =[Wbar;X(:,:,i)-mean(X(:,:,i))];
end
[U,S,V] = svd(Wbar,'econ');
qtilde  = V;
for  j=1:pmax
    x = qtilde(:,j);
    qsum =0;
    for r = 1:j-1
        qsum = qsum +(Q(:,r)*(Q(:,r)'*x));
    end
    x = x-qsum;
    for i=1:k
        b=Wi(:,:,i)*x;
        mu(i) = b'*b;
    end 
    x_old =x;
    for i=1:lmax
        W = [];
        for i = 1:k
            W =[W;sqrt(n(i))*Wi(:,:,i)/sqrt(mu(i))];
        end
        a    = W*x;
        b    = W'*a;
        qsum = 0;
        for r = 1:j-1
        qsum = qsum +(Q(:,r)*(Q(:,r)'*b));
        end
        y=b - qsum;
        x=y/sqrt(y'*y);
        for i=1:k
            b=Wi(:,:,i)*x;
            mu(i) = b'*b;
        end
    xerr = (x-x_old)./x_old;
    mxerr = max(abs(xerr(:)));
    x_old = x;
    end%lmax
    Q(:,j)=x;
end %pmax
for j =1:pmax
    for i=1:k
        c=Wi(:,:,i)*Q(:,j);
        Lambda(j,i) = c'*c;
    end
end
end