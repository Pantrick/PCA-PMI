function cmiv=pmigauss(v1,v2,vcs)
    v1=reshape(v1,[],1);
    v2=reshape(v2,[],1);
 if  nargin==2
        c1=det(cov(v1));
        c2=det(cov(v2));
        c3=det(cov(v1,v2));
        cmiv=0.5*log(c1*c2/c3); 
     elseif  nargin==3
       n1 = size(vcs,2);
       n = n1 +2;

       Cov1 = var(v1);
       Cov2 = var(v2);
       Covm = cov([v1,v2,vcs]);
       Covm1 = cov(vcs);
       Covm2 =cov([v1,vcs]);
       Covm3 = cov([v2,vcs]);

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
       cmiv = 0.5 * (trace(InvC*Covm)+log(C0)-n); 
       
 end
    % cmiv=abs(cmiv);
      if  cmiv==inf 
            cmiv=0;
     end
end
