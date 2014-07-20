function savefig(fig,fname)
% save figures with big fonts and thick lines
%
% use:
%   savefig(figs,fname)
%
% input:
%  figs  - array of figure numbers
%  fname - filename, figures will be saved as filename_{a,b,c,...}
%

%% parse filename
tmp=strfind(fname,'/');
if isempty(tmp)
label = fname;
else
label = fname(tmp(end)+1:end);
end

n = length(fig);

%% save figures
for k = 1:n
    figure(fig(k));
    set(get(gca,'Xlabel'),'fontsize',20);
    set(get(gca,'Ylabel'),'fontsize',20);
    set(gca,'fontsize',20);
    set(get(gca,'title'),'string',[]);
    try 
        set(get(gca,'Children'),'linewidth',2);
    catch
       % 
    end
    print(fig(k),'-depsc',[fname '_' num2str(char(96+k))]);
    close(fig(k));
end

