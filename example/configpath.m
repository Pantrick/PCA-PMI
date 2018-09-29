curpath = pwd;

if curpath(end) == filesep && curpath == ''
else
   curpath = [curpath,filesep];
end
DATA10=[curpath,'data10/'];
DATA50=[curpath,'data50/'];
DATA100=[curpath,'data100/'];
GOLD=[curpath,'gold/'];
RESULT=[curpath,'result/'];
LIB=[curpath,'../lib'];
CODE=curpath;
% if isempty(dir('path.mat'))
%     DATA10=input('Please input the directory of data10(with quotation '' in the begining and end)\n:');
%     DATA50=input('Please input the directory of data50(with quotation '' in the begining and end)\n:');
%     DATA100=input('Please input the directory of data100(with quotation '' in the begining and end)\n:');
%     GOLD=input('Please input the directory of gold(with quotation '' in the begining and end)\n:');
%     CODE=input('Please input the directory of codes(with quotation '' in the begining and end)\n:');
% end
if isunix
    for i={DATA10,DATA50,DATA100,GOLD}
        cd(i{:});
        [~,~]=system('rename .tsv .txt *.tsv');
    end
end
cd(CODE)