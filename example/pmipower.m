%Test the power of PMI in different signal noise ratial
%2015.10.09
SAVEFIG = 0;
%%
T=1000;
n=1000;
Repeat=5;
Amplitude=logspace(log(2^(-2))/log(10),log(8)/log(10),11);
SR=0.9;%Signal Ratio from Z to X, you can choose [0 0.85 0.95 1]
Sxz=0.5;%Signal Ratio from X Z to Y

%%
N = floor((n.*exp(1))^(1/4));

fName={'Linear';'Quadratic';'Cubic';'Sinusoidal';'Exponential';'Checkerboard';'Circular';'CrossShaped';'Sigmoid';'Random'};
fx={@(x)x,...
    @(x)x.^2,...
    @(x)x.*(x-2).*(x+2),...
    @(x)sin(pi.*x),...
    @(x)2.^x,...
    @(x)(mod(x,2)>1)+randi([floor(min(reshape(x,[],1))/2),floor(max(reshape(x,[],1))/2)],size(x))*2+rand(size(x)),...
    @(x)(-1).^binornd(1,0.5,size(x)).*sqrt(max(abs(reshape(x,[],1))).^2-x.^2),...
    @(x)(-1).^binornd(1,0.5,size(x)).*x,...
    @(x)exp(pi.*x)./(1+exp(pi.*x)),...
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
%%
POWERpmi=nan(numel(fx),numel(Amplitude));
% POWERqmi=nan(numel(fx),numel(Amplitude));
POWERcmi=nan(numel(fx),numel(Amplitude));
for i=1:numel(fx)
    parfor j=1:numel(Amplitude)
        
        Z = random('unif',-1,1,n,T);
        Residual = random('unif',-1,1,n,T);
        %X=(1-SR).*Residual+SR.*Z;
        %X=Residual+0.*Z;
        X=(0.5).*Residual+(0.5).*Z;
        %X=0.1*Residual+0.9*Z;
        yeta = randn(n,T);
        
        % Y=SX+SZ=f(X)+f(Z)
        % Ypsu=f(Xpsu)+f(Z);
        SX=feval(fx{i},X);%Signal=X(x)

        SZ=feval(fz{i},Z);%Signal=Z(z)
        %Yst=(SX+SZ)*diag(1./std(SX+SZ))+(yeta)*diag(1./std(yeta))*Amplitude(j);%Y=Signal+A*yeta
        Yst=(SX.*Sxz+SZ.*(1-Sxz))*diag(1./std(SX.*Sxz+SZ.*(1-Sxz)))+(yeta)*diag(1./std(yeta))*Amplitude(j);
        Yrand=cell2mat(arrayfun(@(I,J)Yst(randperm(size(Yst,1))',I),repmat(1:T,[1,1,Repeat]),repmat(reshape(1:Repeat,1,1,[]),[1,T,1]),'UniformOutput',0));
        Xrand=cell2mat(arrayfun(@(I,J)X(randperm(size(X,1))',I),repmat(1:T,[1,1,Repeat]),repmat(reshape(1:Repeat,1,1,[]),[1,T,1]),'UniformOutput',0));
        % Zrand=cell2mat(arrayfun(@(I,J)Z(randperm(size(Z,1))',I),repmat(1:T,[1,1,Repeat]),repmat(reshape(1:Repeat,1,1,[]),[1,T,1]),'UniformOutput',0));
        
        PMIh1=arrayfun(@(I)(pmi(X(:,I),Yst(:,I),Z(:,I),N)),1:T);
        PMIh0_y=arrayfun(@(I,J)pmi(X(:,I),Yrand(:,I,J),Z(:,I),N),repmat(1:T,[1,1,Repeat]),repmat(reshape(1:Repeat,1,1,[]),[1,T,1]));
        PMIh0_x=arrayfun(@(I,J)pmi(Xrand(:,I,J),Yst(:,I),Z(:,I),N),repmat(1:T,[1,1,Repeat]),repmat(reshape(1:Repeat,1,1,[]),[1,T,1]));
        
%         QMIh1=arrayfun(@(I)(qmi(X(:,I),Yst(:,I),Z(:,I),N)),1:T);
%         QMIh0_y=arrayfun(@(I,J)qmi(X(:,I),Yrand(:,I,J),Z(:,I),N),repmat(1:T,[1,1,Repeat]),repmat(reshape(1:Repeat,1,1,[]),[1,T,1]));
%         QMIh0_x=arrayfun(@(I,J)qmi(Xrand(:,I,J),Yst(:,I),Z(:,I),N),repmat(1:T,[1,1,Repeat]),repmat(reshape(1:Repeat,1,1,[]),[1,T,1]));
%         
%         
        CMIh1=arrayfun(@(I)(cmi(X(:,I),Yst(:,I),Z(:,I),N)),1:T);
        CMIh0_y=arrayfun(@(I,J)cmi(X(:,I),Yrand(:,I,J),Z(:,I),N),repmat(1:T,[1,1,Repeat]),repmat(reshape(1:Repeat,1,1,[]),[1,T,1]));
        CMIh0_x=arrayfun(@(I,J)cmi(Xrand(:,I,J),Yst(:,I),Z(:,I),N),repmat(1:T,[1,1,Repeat]),repmat(reshape(1:Repeat,1,1,[]),[1,T,1]));
        
        PMIh0=[PMIh0_y(:);PMIh0_x(:)];
        PMIh0=sort(reshape(PMIh0,[],1),'descend');
        ALPHAp=mean(PMIh0([ceil(numel(PMIh0)*0.05) floor(numel(PMIh0)*0.05)]));
        POWERpmi(i,j)=sum(PMIh1>=ALPHAp)/numel(PMIh1);
        
%         QMIh0=[QMIh0_y(:);QMIh0_x(:)];
%         QMIh0=sort(reshape(QMIh0,[],1),'descend');
%         ALPHAq=mean(QMIh0([ceil(numel(QMIh0)*0.05) floor(numel(QMIh0)*0.05)]));
%         POWERqmi(i,j)=sum(QMIh1>=ALPHAq)/numel(QMIh1);
        
        CMIh0=[CMIh0_y(:);CMIh0_x(:)];
        CMIh0=sort(reshape(CMIh0,[],1),'descend');
        ALPHAc=mean(CMIh0([ceil(numel(CMIh0)*0.05) floor(numel(CMIh0)*0.05)]));
        POWERcmi(i,j)=sum(CMIh1>=ALPHAc)/numel(CMIh1);
        [i j]
    end
end

%%
figure(1);
imagesc(POWERcmi);
%Position=[253 41 849 728];
%set(gcf,'Outerposition',Position);
%fName={'Linear';'Quadratic';'Sinusoid';'Exponential';'Cubic';'Random'};
set(gca,'YTick',1:10,'YTickLabel',fName);
set(gca,'FontSize',15,'FontWeight','bold');
axis square;
axis tight;

colormap(redgreencmap(256, 'Interpolation','Linear'));
colorbar('FontSize',10,'FontWeight','bold');
xlabel('Noise Amplitude','FontSize',20,'FontWeight','bold');
set(gca,'XTickLabel',{0.5,1,2,4,8},'Clim',[0,1]);
title({'Power of CMI'; '_{with noises in different Amplitude}'},'FontSize',18,'FontWeight','bold');


figure(2)
imagesc(POWERpmi);
%Position=[253 41 849 728];
%set(gcf,'Outerposition',Position);
%fName={'Linear';'Quadratic';'Sinusoid';'Exponential';'Cubic';'Random'};
set(gca,'YTick',1:10,'YTickLabel',fName);
set(gca,'FontSize',15,'FontWeight','bold');
axis square;
axis tight;

colormap(redgreencmap(256, 'Interpolation','Linear'));
colorbar('FontSize',10,'FontWeight','bold');
xlabel('Noise Amplitude','FontSize',20,'FontWeight','bold');
set(gca,'XTickLabel',{0.5,1,2,4,8},'Clim',[0,1]);
title({'Power of PMI'; '_{with noises in different Amplitude}'},'FontSize',18,'FontWeight','bold');

figure(3);
imagesc([POWERpmi;POWERcmi])
%fName={'Linear';'Quadratic';'Sinusoid';'Exponential';'Cubic';'Random'};
set(gca,'YTick',1:20,'YTickLabel',{fName{:},fName{:}})
set(gca,'FontSize',15,'FontWeight','bold');
axis square;
axis tight;

colormap(redgreencmap(256, 'Interpolation','Linear'));
colorbar('FontSize',10,'FontWeight','bold');
xlabel('Noise Amplitude','FontSize',20,'FontWeight','bold');
set(gca,'XTickLabel',{0.5,1,2,4,8},'Clim',[0,1]);
ylabel('CMI          PMI','FontSize',40,'FontWeight','bold');

% figure();
% imagesc(POWERpdcor);
% %Position=[253 41 849 728];
% %set(gcf,'Outerposition',Position);
% %fName={'Linear';'Quadratic';'Sinusoid';'Exponential';'Cubic';'Random'};
% set(gca,'YTick',1:10,'YTickLabel',fName);
% set(gca,'FontSize',15,'FontWeight','bold');
% axis square;
% axis tight;

colormap(redgreencmap(256, 'Interpolation','Linear'));
colorbar('FontSize',10,'FontWeight','bold');
xlabel('Noise Amplitude','FontSize',20,'FontWeight','bold');
set(gca,'XTickLabel',{0.5,1,2,4,8},'Clim',[0,1]);
title({'Power of Pdcor'; '_{with noises in different Amplitude}'},'FontSize',18,'FontWeight','bold');

%%
if SAVEFIG == 1
    cd(RESULT)
    saveas(figure(3),'power_depend.png','png')
    saveas(figure(3),'power_depend.eps','psc2')
%     saveas(figure(3),'power_depend.emf','emf')
%     saveas(figure(3),'power_depend.tif','tif')
    saveas(figure(1),'powercmi_depend.png','png')
    saveas(figure(1),'powercmi_depend.eps','psc2')
%     saveas(figure(1),'powercmi_depend.emf','emf')
%     saveas(figure(1),'powercmi_depend.tif','tif')
    saveas(figure(2),'powerpmi_depend.png','png')
    saveas(figure(2),'powerpmi_depend.eps','psc2')
%     saveas(figure(2),'powerpmi_depend.emf','emf')
%     saveas(figure(2),'powerpmi_depend.tif','tif')
end