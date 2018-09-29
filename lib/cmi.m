function [ I E ] = cmi( x,y,z,n )
%2nd edition start on 2014/5/13
%using F instead of P , P=F/N;
%Conditional Mutual Information= sum( P(x,y,z)*log(P(x,y,z)*P(z)/(P(x,z)*P(y,z)))  );
%also = sum( F(x,y,z)/N*log(F(x,y,z)*F(z)/(F(x,z)*F(y,z)))  );
if nargin==3
    n=10;
end
if numel(n)==1
    Fxyz=hist4([x,y,z],[n,n,n]);
elseif numel==3
    Fxyz=hist4([x,y,z],n);
else
    error(message('stats:cmi:Number of N must be 1 or 3.'));
end

%Pxyz=Fxyz./size(x,1);
N=size(x,1);  
Fxz=sum(Fxyz,2);
Fyz=sum(Fxyz,1);
Fz=sum(sum(Fxyz));
Fz_FxzFyz=cellfun(@(z,x,y)(z./(x*y)),num2cell(Fz,[1 2]),num2cell(Fxz,1),num2cell(Fyz,2),'UniformOutput',0);
Fz_FxzFyz=cell2mat(Fz_FxzFyz);
E=Fxyz./N.*log(Fxyz.*Fz_FxzFyz);
E(Fxyz==0)=0;
I=sum(sum(sum(E)));
end

