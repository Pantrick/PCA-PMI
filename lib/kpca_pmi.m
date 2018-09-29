function [G,Gval,order]=kpca_pmi(data,lambda,order0)
% KPCA_PMI is an function to construct the Gene Expression Network by using
% Pathing-consist Algorithm with PMI in Kernel Version.
% The "Distance Correlation" is used as the Kernel in this function.
% G=KPCA_PMI(data,lambda) "data" is a n*T matrix, with n genes and T samples,
% the network of the n genes will be constructed. "lambda" is a threshold
% of reducing edges.
% G=KPCA_PMI(data,lambda,order0) "order0" is a positive integer, if the
% algorithm is processing in the order larger than "order0", the algorithm
% will stop reducing the edge and return the network.
% G=KPCA_PMI(...) "G" returns a binary lower triangle square matrix, the
% element "1" means there is an edge in the matching place, while the
% element "0" means there is no edge in the matching place in the network.
% [G,GVAL]=KPCA_PMI(...) "Gval" returns a lower triangle square matrix,
% every entry in "Gval" represents the strength of the matching edge
% [G,GVAL,order]=KPCA_PMI(...) "order" returns a final order where the
% algorithm processes

    n_gene=size(data,1); %gene number
    G=ones(n_gene,n_gene);
    G=tril(G,-1)'; 
    G=G+G';
    Gval=G;
    order=-1;t=0;
    if nargin==2 || isempty(order)
        order0=n_gene;
    end
    corrmat=dcorr(data'); % compute the Kernel Matrix

    while t==0
        order=order+1;
        if order>order0
           G=tril(G,-1)';     Gval=tril(Gval,-1)';  
           order=order-1; % The value of order is the last order of pc algorith 
           return
        end

        [G,Gval,t]=edgereduce(G,Gval,order,corrmat,t,lambda); %reduce edge

         if t==0
              %disp('No edge is reduce! Algorithm  finished!'); %For Debugging
              break;
         else 
              t=0;
         end
         %disp(t);%For Debugging
    end
  
   G=tril(G,-1)';     Gval=tril(Gval,-1)';  
   order=order-1; % The value of order is the last order of pc algorith 
   
end

%% edgereduce is pca_pmi
function [G,Gval,t]=edgereduce(G,Gval,order,corrmat,t,lambda)
    Gprior=G;
    if order==0
        Gval=-1/2.*log(1-(corrmat-eye(size(corrmat,1))).^2);
        G(Gval<lambda)=0;
        t=t+1;
    else

        [x,y]=find(triu(G,1)~=0);
        for i=1:numel(x)
    %               i;%For debugging
    %               j;%For debugging
            adj=find(Gprior(x(i),:)&Gprior(y(i),:)); 
            if size(adj,2)>=order
                %combntnslist=combntns(adj,order); %For Debugging
                combntnslist=nchoosek(adj,order);
                combntnsrow=size(combntnslist,1);   
                pmiv=arrayfun(@(k)pmi(corrmat,x(i),y(i),reshape(combntnslist(k,:),1,[])),1:combntnsrow);
                %pmiv(isinf(pmiv))=0;
                pmiv=min(pmiv);
    %                    pmiv;%For debugging
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
       % C = inv(InvC);%For debugging  

       C0 = (det(Covm)*det(Covm1)/(det(Covm2)*det(Covm3)))*Cov1*Cov2*(InvCovm(2,2)-InvCovm3(1,1)+InvCov2)*(InvCovm(1,1)-InvCovm2(1,1)+InvCov1);
       pmiv = 0.5 * (trace(InvC*Covm)+log(C0)-n); 
       
 end
     if  pmiv==inf 
            pmiv=0;
     end
       
end


