% This is matlab code for conditional mutual information based path
% consitency algorithm for inferring network by the conditional
% dependence. 
% Input:
% 'data' is expression of variable,in which row is varible and column is the sample;
% 'lamda' is the parameter decide the dependence;
% 'order0' is the parameter to end the program when order=order0;
% Output:
% 'G' is the 0-1 network or graph after pc algorith
% 'Gval' is the network with strenthness of dependence;
% 'order' is the order of the pc algorithm, here is equal to order0;
% If nargin==2,the algorithm will be terminated untill there is no higher
% order networks.
% This version is revised for search the circle in the network G;
% Version Data: 2011-7-7

function [G,Gval,order]=pca_cmi(data,lamda,order0)  
                                               
n_gene=size(data,1);
G=ones(n_gene,n_gene);
G=tril(G,-1)'; 
G=G+G';
Gval=G;
order=-1;t=0;
while t==0
     order=order+1;
     if nargin==3
       if order>order0
           G=tril(G,-1)';     Gval=tril(Gval,-1)';  
           order=order-1; % The value of order is the last order of pc algorith 
           return
       end
     end
    [G,Gval,t]=edgereduce(G,Gval,order,data,t,lamda);
 
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

%% edgereduce is pca_cmi
function [G,Gval,t]=edgereduce(G,Gval,order,data,t,lamda)
%[nrow,ncol]=find(G~=0);

if order==0
    for i=1:size(G,1)-1
        for j=i+1:size(G,1)
            if G(i,j)~=0
                cmiv=cmi(data(i,:),data(j,:));
                Gval(i,j)=cmiv;  Gval(j,i)=cmiv;
                if cmiv<lamda
                    G(i,j)=0;G(j,i)=0;
                end
            end
        end
    end
          t=t+1;
else
  for i=1:size(G,1)-1
      for j=i+1:size(G,1)
          if G(i,j)~=0
%               i ;
%               j;
              adj=[] ;
              for k=1:size(G,1)
                  if G(i,k)~=0 && G(j,k)~=0
                      adj=[adj,k] ;
                  end
              end
              if size(adj,2)>=order
                   combntnslist=combntns(adj,order);
                   combntnsrow=size(combntnslist,1);   
                   cmiv=0;
                    v1=data(i,:);v2=data(j,:);
                   for k=1:combntnsrow   
                    vcs=data(combntnslist(k,:),:);  %vcs=data(adj(k),:);
                     a=cmi(v1,v2,vcs) ;
                     cmiv=max(cmiv,a);
                   end
%                    cmiv;
                      Gval(i,j)=cmiv; Gval(j,i)=cmiv;
                      if cmiv<lamda
                          G(i,j)=0; G(j,i)=0;
                      end              
                           t=t+1;
              end
          end
                      
      end
  end
  
   
               
end
end

%% compute conditional mutual information of x and y 
function cmiv=cmi(v1,v2,vcs)
 if  nargin==2
        c1=det(cov(v1));
        c2=det(cov(v2));
        c3=det(cov(v1,v2));
        cmiv=0.5*log(c1*c2/c3); 
     elseif  nargin==3
        c1=det(cov([v1;vcs]'));
        c2=det(cov([v2;vcs]'));
        c3=det(cov(vcs'));
        c4=det(cov([v1;v2;vcs]'));
        cmiv=0.5*log((c1*c2)/(c3*c4));       
 end
    % cmiv=abs(cmiv);
     if  cmiv==inf 
            cmiv=0;
     end

end
