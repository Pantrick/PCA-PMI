%pcalg_partcor is path consistency algorithm based on partial correlation
%input variable is expression data and delete edge parameter:data and lamda
%first generate a whole graph G incording to the size of genes from data
%sencond compute the partial correlation throgh function
%partial_correlation(var_x,var_y,var_z)
%thirdly reduce genes based on the the partial correlation\
% example:
% lamda=0.1;
% [G,Gval,l]=pcalg_partcor(data,lamda);
%% dream3 data/ compute the tpr,fpr,ppv,precesion
 

function [G,Gval,order]=pcapcc(data,lamda,order0)
n_gene=size(data,1);
G=ones(n_gene,n_gene);
G=tril(G,-1)';  Gval=zeros(size(G));
C=corrcoef(data');           % C is the correalton coefficiet matrix
% l=-1;t=0;
% while t==0
%     l=l+1 
%     [G,Gval,t]=edgereduce_pcc(G,Gval,l,C,t,lamda) 
%      if t==0
%           disp('No edge is reduce! Algorithm  finished!');
%           break;
%      else 
%           t=0;
%      end
% end
% 
% end
order=-1;t=0;
while t==0
     order=order+1;
     if nargin==3
       if order>order0
%            G=tril(G,-1)';     Gval=tril(Gval,-1)';  
           order=order-1; % The value of order is the last order of pc algorith 
           return
       end
     end
    [G,Gval,t]=edgereduce_pcc(G,Gval,order,C,t,lamda) ;
 
     if t==0
          disp('No edge is reduce! Algorithm  finished!');
          break;
     else 
          t=0;
     end
end

%% edgereduce is reduce edge incording to the partial corrlation coeffie
function [G,Gval,t]=edgereduce_pcc(G,Gval,order,C,t,lamda)
[nrow,ncol]=find(G~=0);
if order==0
  for m=1:size(nrow,1)
          i=nrow(m);j=ncol(m);
          %v1=data(i,:);v2=data(j,:);
          partcorcoef=partial_correlation_coefficient(C,i,j);
              Gval(i,j)=partcorcoef;         
          if abs(partcorcoef)<lamda
            G(i,j)=0;
          
          end
  end
          t=t+1;
else
  for m=1:size(nrow,1)
        i=nrow(m);j=ncol(m);
        adj=find(G(i,:)~=0);
        adj(adj<=j)=[];
        size_adj=size(adj,2);
          if size_adj>=order                
            combntnslist=combntns(adj,order);
            combntnsrow=size(combntnslist,1);
            combntnslist=combntnslist';
            pcc=0;
               for k=1:combntnsrow   
                 % v1=data(i,:);v2=data(j,:);vcs=data(combntnslist(k,:),:);  %vcs=data(adj(k),:);
                  %cmiv=cmi(v1,v2,vcs);
                   
                   partcorcoef=partial_correlation_coefficient(C,i,j,combntnslist(:,k));
                  
                   %Gval(i,j)=partcorcoef;  
                   if abs(partcorcoef)<lamda
                      G(i,j)=0;
                      Gval(i,j)=partcorcoef;
                     % S(i,j)=adj(k);
                   end           
               end
               t=t+1;
           end           
          
   end 
end

% function [G,Gval,t]=edgereduce_pcc(G,Gval,l,C,t,lamda)
% [nrow,ncol]=find(G~=0);
% if l==0
%   for m=1:size(nrow,1)
%           i=nrow(m);j=ncol(m);
%           %v1=data(i,:);v2=data(j,:);
%           partcorcoef=partial_correlation_coefficient(C,i,j);
%               Gval(i,j)=partcorcoef;         
%           if abs(partcorcoef)<lamda
%             G(i,j)=0;
%           
%           end
%   end
%           t=t+1;
% else
%   for m=1:size(nrow,1)
%         i=nrow(m);j=ncol(m);
%         adj=find(G(i,:)~=0);
%         adj(adj<=j)=[];
%         size_adj=size(adj,2);
%           if size_adj>l                
%             combntnslist=combntns(adj,l);
%             combntnsrow=size(combntnslist,1);
%             combntnslist=combntnslist';
%                for k=1:combntnsrow   
%                  % v1=data(i,:);v2=data(j,:);vcs=data(combntnslist(k,:),:);  %vcs=data(adj(k),:);
%                   %cmiv=cmi(v1,v2,vcs);
%                    
%                    partcorcoef=partial_correlation_coefficient(C,i,j,combntnslist(:,k));
%                   
%                    %Gval(i,j)=partcorcoef;  
%                    if abs(partcorcoef)<lamda
%                    G(i,j)=0;
%                  Gval(i,j)=partcorcoef;
%                      % S(i,j)=adj(k);
%                    end           
%                end
%                t=t+1;
%            end           
%           
%    end 
% end
               
end
 
end
%% partial correlation computer the conditional correlation between 
%varible x and variable y under condition variables z

function [partcorcoef]=partial_correlation_coefficient(C,i,j,condition_k)
 
if nargin==3
        partcorcoef=C(i,j);
end
if nargin==4
    n=size(condition_k,2);
    if n==0
        partcorcoef=partial_correlation_coefficient(C,i,j);
    else
        k=condition_k(n);
        condition_k(n)=[];
        
        a=partial_correlation_coefficient(C,i,j,condition_k);
        b=partial_correlation_coefficient(C,i,k,condition_k);
        c=partial_correlation_coefficient(C,j,k,condition_k);
        partcorcoef=(a-b*c)/sqrt((1-b^2)*(1-c^2));
%          partcorcoef=abs(partcorcoef);
    end
    
end


end
