% % We consider the following relationships between X and Y when given Z for
% % simulation study
% %'Linear';'Quadratic';'Cubic';'Sinusoidal';'Exponential';'Checkerboard';'Circular';'CrossShaped';'Sigmoid';'Random'
% % Showing mode can be chose from 'clean' or 'noise'
%%
n = 5000;
SAVEFIG = 0;
MODE = 'noise'; %MODE can be chose from 'clean' or 'noise'
%%
fName={'Linear';'Quadratic';'Cubic';'Sinusoidal';'Exponential';'Checkerboard';'Circular';'CrossShaped';'Sigmoid';'Random'};
Xx={@(n)linspace(-3,3,n)';
    @(n)linspace(-3,3,n)';
    @(n)linspace(-3,3,n)';
    @(n)linspace(-3,3,n)';
    @(n)linspace(-3,3,n)';
    @(n)linspace(-3,3,n)';
    @(n)linspace(-3,3,n)';
    @(n)linspace(-3,3,n)';
    @(n)linspace(-3,3,n)';
    @(n)random('normal',0,1,[1,n])};...
    
fx={@(x)x,
    @(x)x.^2,
    @(x)x.*(x-2).*(x+2),
    @(x)sin(pi.*x),
    @(x)2.^x,
    @(x)(mod(x,2)>1)+randi([floor(min(reshape(x,[],1))/2),floor(max(reshape(x,[],1))/2)],size(x))*2+rand(size(x)),
    @(x)(-1).^binornd(1,0.5,size(x)).*sqrt(max(abs(reshape(x,[],1))).^2-x.^2),
    @(x)(-1).^binornd(1,0.5,size(x)).*x,
    @(x)exp(pi.*x)./(1+exp(pi.*x))
    @(x)random('normal',0,1,size(x))};

Noisex={@(x)random('normal',0,std(x)/5,size(x));
    @(x)random('normal',0,std(x)/5,size(x));
    @(x)random('normal',0,std(x)/5,size(x));
    @(x)random('normal',0,std(x)/5,size(x));
    @(x)random('normal',0,std(x)/5,size(x));
    @(x)0;
    @(x)random('normal',0,std(x)/5,size(x));
    @(x)random('normal',0,std(x)/5,size(x));
    @(x)random('normal',0,std(x)/5,size(x));
    @(x)0};...
position=[1,2,3,4,5,6,7,8,9,11];
x=cellfun(@(f)feval(f,n),Xx,'UniformOutput',0);
y=cellfun(@(f,X)feval(f,X),fx,x,'UniformOutput',0);
if ~strcmp(MODE,'clean') && strcmp(MODE,'noise')
    y=cellfun(@(Noise,X)feval(Noise,X)+X,Noisex,y,'UniformOutput',0);
end
for i=1:numel(fx) ...
    subplot(4,3,position(i));
    plot(x{i},y{i},'.','MarkerSize',0.5);
    hold on;
    axis([-3,3,-3,3]);
    axis tight;
    %axis equal;
    axis off;...
    %box off;
    title(fName{i});
end
if SAVEFIG == 1
    cd(RESULT)
    saveas(gcf,sprintf('funxy_%s_axisoff.png',MODE),'png');
    saveas(gcf,sprintf('funxy_%s_axisoff.png',MODE),'psc2');
end