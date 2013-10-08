% function plot_abundance(varargin)
function plot_abundance(f,tit,outname)

%    varargin: list of input files to plot on abundance plot

%c={'bo','ro','go','co','ko','ks','bs','rs','gs','cs'};
c={'bo','ro','go','co','ko'};
%c={'ks','bs','rs','gs','cs'};

%leg={'original','post correction'};
leg={'original','matches only','post correction'};

fig=figure;
for i=1:length(f)
    [pathstr,name,ext] = fileparts(f{i});
    fin=fopen(f{i},'r');
    data=textscan(fin,'%d%d');
    nunq=data{1};
    uc=data{2};
    fclose(fin);
    
    loglog(uc,nunq,c{i},'MarkerSize',6)
    hold on
end
    title(tit,'FontSize',16)
    xlabel('copy number (reads)','FontSize',16)
    ylabel('unique CDR3/VJ combination','FontSize',16)
    set(gca,'FontSize',14) 
    legend(leg,'FontSize',14)

    %outdir='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/analysis/plots/abundance/';
    %outpath=[outdir,outname];
    print(fig,'-dpdf','-r300',[outname,'.pdf']);

