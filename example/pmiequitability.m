%20150917
% Test the equitability of PMI
%% Parameter Setting
n=500;%sample size
T=100;%repeat for each R2
SAVEFIG = 0;

%%
rng(cputime*10);
R2=[1:-0.01:0] ;
R_square = zeros(6,T,numel(R2));
test_matrix=zeros(6,T,numel(R2));
%Position=[253 41 849 728];

%%
Z = random('Unif',-2,2,n,T) ;
X = random('Unif',-2,2,n,T);
yeta = randn(n,T) ;
yeta=yeta*diag(1./std(yeta));
%%
% fx={@(x)x+1,
% @(x)(x-0.5).*(x-1)+1,
% @(x)sin(x)+pi,
% @(x) 0.5*exp(x)+3,
% @(x) (x+0.5).*(x-.5).*(x-0.25)+0.25,
% @(x)random('Unif',3,4,size(x))+10};
% fz={@(z)z.*z,
% @(z)z.*z,
% @(z)z.*z,
% @(z)z,
% @(z)z.*z,
% @(z)z/20};
% fName={'Linear';'Quadratic';'Sinusoid';'Exponential';'Cubic';'Random'};
fName={'Linear';'Quadratic';'Cubic';'Sinusoidal';'Exponential';'Checkerboard';'Circular';'CrossShaped';'Sigmoid';'Random'};
fx={@(x)x,
    @(x)x.^2,
    @(x)x.*(x-2).*(x+2),
    @(x)sin(pi.*x),
    @(x)2.^x,
    @(x)(mod(x,2)>1)+randi([floor(min(reshape(x,[],1))/2),floor(max(reshape(x,[],1))/2)],size(x))*2+rand(size(x)),
    @(x)(-1).^binornd(1,0.5,size(x)).*sqrt(max(abs(reshape(x,[],1))).^2-x.^2),
    @(x)(-1).^binornd(1,0.5,size(x)).*x,
    @(x)exp(pi.*x)./(1+exp(pi.*x))
    @(x)10.*randn(size(x))};
fz={@(z)z+1;
    @(z)z+1;
    @(z)z+pi;
    @(z)z+3;
    @(z)-z+2.5;
    @(z)0.5*z;
    @(z)z+4;
    @(z)0.5*z;
    @(z)-z;
    @(z)z/20+10};
N  = [repmat(ceil(n^(1/4)),numel(fx),1)] ;
%%
Signal=cell(numel(fx));
for i=1:numel(R2)
    Signal=cellfun(@(f)(feval(f,X)),fx,'UniformOutput',0);%Signal=X(x)
    SZ=cellfun(@(f)(feval(f,Z)),fz,'UniformOutput',0);%Signal=Z(x)
    Y=cellfun(@(s,z)( sqrt(R2(i))*(s+z)+sqrt(1-R2(i))*(yeta)*diag(std(s+z))*diag(1./std(yeta))),Signal,SZ,'UniformOutput',0);%Y=Signal+(1-R2)/R2*yeta
    %for j=3
    parfor j=1:numel(fx)
        test_matrix(j,:,i) = arrayfun(@(I)pmi(X(:,I),Y{j}(:,I),Z(:,I),N(j)),1:T);%pmi(X,Y,Z,N)
        R_square(j,:,i) =arrayfun(@(I)(1-corr(Signal{j}(:,I)+SZ{j}(:,I),Y{j}(:,I)).^2),1:T);
    end
end
%%
COLOR=hsv(numel(fx));
%figure('Outerposition',Position);
figure(1);
set(gcf,'paperunits','centimeters','papersize',[4.1 4.5]);
mean_testmatrix=cellfun(@(x)(mean(x)),num2cell(test_matrix,2));
h=nan(size(fx));
%for j=3
for j=1:numel(fx)
    plot(reshape(R_square(j,:,:),[],1),reshape(test_matrix(j,:,:),[],1),'.','Color',COLOR(j,:),'MarkerSize',0.2);hold on
    h(j)=plot(1-R2,reshape(mean_testmatrix(j,:,:),1,[]),'-','Color',COLOR(j,:),'LineWidth',2);hold on;
end

%for j=3
% for j=1:numel(fx)
%     h(j)=plot(1-R2,reshape(mean_testmatrix(j,:,:),1,[]),'-','Color',COLOR(j,:),'LineWidth',2);hold on;
% end
legend(h,fName,'FontSize',12);
ylim([0,1.4])
%%
%figure('Outerposition',Position);
figure(2);
set(gcf,'paperunits','centimeters','papersize',[4.1 4.5]);
for j=1:numel(fx)
    plot(1-R2,reshape(mean_testmatrix(j,:,:),1,[]),'Color',COLOR(j,:),'LineWidth',2);hold on;
end
legend(fName,'FontSize',12);
ylim([0,1.4])
%%
%figure('Outerposition',Position);
figure(3);
set(gcf,'paperunits','centimeters','papersize',[4.1 4.5]);
for j=1:numel(fx)
    plot(1,2.0,'o','MarkerEdgeColor','none','MarkerFaceColor',COLOR(j,:),'MarkerSize',6);hold on;
end
for j=1:numel(fx)
    plot(reshape(R_square(j,:,:),[],1),reshape(test_matrix(j,:,:),[],1),'o','MarkerEdgeColor','none','MarkerFaceColor',COLOR(j,:),'MarkerSize',1);hold on;
end
legend(fName,'FontSize',12);
%alpha(0.5);
ylim([0,1.4])
%%
for i=1:3
    figure(i);
    set(gca,'FontSize',12,'TickDir','out');
    xlabel('Noise(1-R^2)','Color','black','FontSize',24);
    ylabel('PMI','Color','black','FontSize',24);
    title(sprintf('N=%d',n),'Color','black','FontSize',26,'FontWeight','bold');
end
%%
if SAVEFIG == 1
    cd(RESULT)
    saveas(figure(3),sprintf('n%d_line.png',n),'png');
    saveas(figure(3),sprintf('n%d_line.eps',n),'psc2');
    % 
    saveas(figure(2),sprintf('n%d_plot2.png',n),'png');
    saveas(figure(2),sprintf('n%d_plot2.eps',n),'psc2');

    % 
    saveas(figure(1),sprintf('n%d_plot.png',n),'png');
    saveas(figure(1),sprintf('n%d_plot.eps',n),'psc2');
end
