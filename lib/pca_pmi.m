%2016_3_14_第三次文章修改时,corr+corrmat版本
function [G,Gval,order]=pca_pmi(data,lambda,order0)  
                                               
    n_gene=size(data,1);
    G=ones(n_gene,n_gene);
    G=tril(G,-1)'; 
    G=G+G';
    Gval=G;
    order=-1;t=0;
    if nargin==2 || isempty(order)
        order0=n_gene;
    end
    corrmat=corr(data');

    while t==0
        order=order+1;
        if order>order0
           G=tril(G,-1)';     Gval=tril(Gval,-1)';  
           order=order-1; % The value of order is the last order of pc algorith 
           return
        end

        [G,Gval,t]=edgereduce(G,Gval,order,corrmat,t,lambda);

         if t==0
              %disp('No edge is reduce! Algorithm  finished!');
              break;
         else 
              t=0;
         end
    end
  
   G=tril(G,-1)';     Gval=tril(Gval,-1)';  
   order=order-1; % The value of order is the last order of pc algorith 
end

%% edgereduce is pca_pmi
function [G,Gval,t]=edgereduce(G,Gval,order,corrmat,t,lambda)
    %[nrow,ncol]=find(G~=0);
    %pmi=@(i,j,k)(-1/2*log(1-det(corrmat(i,j)*corrmat(k,k)-corrmat(k,i)*corrmat(j,k))^2/det(corrmat([i,k],[i,k]))/det(corrmat([j,k],[j,k]))));
    pmi=@(i,j,k)(-1/2*log(det(corrmat([i,j,k],[i,j,k]))*det(corrmat(k,k))/det(corrmat([i,k],[i,k]))/det(corrmat([j,k],[j,k]))));
    %mi=@(i,j)(1/2*log(1-corrmat(i,j)^2));
    if order==0
        Gval=-1/2.*log(1-(corrmat-eye(size(corrmat,1))).^2);
        G(Gval<lambda)=0;
        t=t+1;
    else
        [x,y]=find(triu(G,1)~=0);
        for i=1:numel(x)
    %               i ;
    %               j;
            adj=find(G(x(i),:)&G(y(i),:));
            if size(adj,2)>=order
                combntnslist=nchoosek(adj,order);
                combntnsrow=size(combntnslist,1);   
                pmiv=arrayfun(@(k)pmi(x(i),y(i),reshape(combntnslist(k,:),1,[])),1:combntnsrow);
                pmiv(isinf(pmiv))=0;
                pmiv=max(pmiv);
    %                    pmiv;
                Gval(x(i),y(i))=pmiv; Gval(y(i),x(i))=pmiv;
                if pmiv<lambda
                  G(x(i),y(i))=0; G(y(i),x(i))=0;
                end              
                t=t+1;
            end

        end
    end
end

%% compute partial mutual information of x and y
function pmiv=pmi(corrmat,x,y,z)
 if  nargin==3
        c1=det(corrmat(x,x));
        c2=det(corrmat(y,y));
        c3=det(corrmat([x,y],[x,y]));
        pmiv=0.5*log(c1*c2/c3); 
     elseif  nargin==4
       n1 = numel(z);
       n  = n1 +2;
       z=reshape(z,1,[]);
       Cov1 = corrmat(x,x);
       Cov2 = corrmat(y,y);
       Covm = corrmat([x,y,z],[x,y,z]);
       Covm1 = corrmat(z,z);
       Covm2 =corrmat([x z],[x z]);
       Covm3 = corrmat([y z],[y z]);

       InvCov1 = 1/Cov1;
       InvCov2 = 1/Cov2;
       InvCovm = inv(Covm);
       InvCovm1 = inv(Covm1);
       InvCovm2 = inv(Covm2);
       InvCovm3 = inv(Covm3);

       C11 = -InvCovm(1,2)*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*InvCovm(1,2)+InvCovm(1,1);
       C12 = 0;
       C13 = -InvCovm(1,2)*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*(InvCovm(2,3:2+n1)-InvCovm3(1,2:1+n1))+InvCovm(1,3:2+n1);
       C23 = -InvCovm(1,2)*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*(InvCovm(1,3:2+n1)-InvCovm2(1,2:1+n1))+InvCovm(2,3:2+n1);
       C22 = -InvCovm(1,2)*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*InvCovm(1,2)+InvCovm(2,2);
       C33 = -(InvCovm(2,3:2+n1)-InvCovm3(1,2:1+n1))'*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*(InvCovm(2,3:2+n1)-InvCovm3(1,2:1+n1))-(InvCovm(1,3:2+n1)-InvCovm2(1,2:1+n1))'*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*(InvCovm(1,3:2+n1)-InvCovm2(1,2:1+n1))+(InvCovm(3:2+n1,3:2+n1)-InvCovm3(2:1+n1,2:1+n1))+(InvCovm(3:2+n1,3:2+n1)-InvCovm2(2:1+n1,2:1+n1))+InvCovm1;
       InvC = [[C11,C12,C13];[C12,C22,C23];[[C13',C23'],C33]];
       % C = inv(InvC);  

       C0 = (det(Covm)*det(Covm1)/(det(Covm2)*det(Covm3)))*Cov1*Cov2*(InvCovm(2,2)-InvCovm3(1,1)+InvCov2)*(InvCovm(1,1)-InvCovm2(1,1)+InvCov1);
       pmiv = 0.5 * (trace(InvC*Covm)+log(C0)-n); 
       
 end
     if  pmiv==inf 
            pmiv=0;
     end
       
end

