function PDCOR=pdcor(X,Y,Z)
% 2nd version, can be used by GPU
% PDCOR=pdcor(X,Y,Z)
% pdcor() calculates the Partial Distance Correlation Between X and Y conditioned on Z
% Requirement: Matlab Statistics Toolbox
% X, Y and Z should have the same number of rows.
% If Z is emtpy or there are only two input arguments
% pdcor() calculates the dcor() between X and Y
%reference ____
if size(X,1)~=size(Y,1) | size(X,1)~=size(Z,1)
    error('X, Y and Z should have the same number of rows.');
end
    N=size(X,1);
    Xij=samplenorm(X,N);
    Yij=samplenorm(Y,N);
    Zij=samplenorm(Z,N);
    %Xij=squareform(tril(Xij,-1));
    %Yij=squareform(tril(Yij,-1));
    %Zij=squareform(tril(Zij,-1));
    normZij=sum(Zij.*Zij);
    X_z=Xij-sum(Xij.*Zij)*Zij./normZij;
    Y_z=Yij-sum(Yij.*Zij)*Zij./normZij;
    PDCOR=sum(X_z.*Y_z)./sqrt(sum(X_z.^2))./sqrt(sum(Y_z.^2));
    PDCOR=gather(PDCOR);
    function Aij=samplenorm(data,N)
        %aij=squareform(pdist(data));
        Xz=reshape(data,N,1,[]);
        XzT=reshape(data,1,N,[]);
        aij=sqrt(sum(bsxfun(@minus,Xz,XzT).^2,3));
        Aij=aij-ones(N,1)*sum(aij,1)/(N-2)-sum(aij,2)*ones(1,size(aij,1))/(N-2)+sum(aij(:))/(N-1)/(N-2);
        %Aij=reshape(Aij,[],1);
        Aij=Aij(tril(true(N),-1));
        Aij=Aij(:);
    end
end




