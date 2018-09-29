function corrmat=dcorr(data)
%DCOV Summary of this function goes here
%   Detailed explanation goes here
    ddata=arrayfun(@(i)samplenorm(data(:,i)),1:size(data,2),'UniformOutput',0);
    ddata=cell2mat(ddata);
    N=size(data,1);
    covmat=sqrt(ddata'*ddata/(N^2));
    corrmat=sqrt(diag(1./diag(covmat)))*covmat*sqrt(diag(1./diag(covmat)));
    function Aij=samplenorm(data)
        data=reshape(data,[],1);
        aij=abs(data*ones(1,size(data,1))-ones(size(data,1),1)*data');
        Aij=aij-ones(size(aij,1),1)*mean(aij,1)-mean(aij,2)*ones(1,size(aij,1))+mean(mean(aij));
        Aij=reshape(Aij,[],1);
    end
end


