% function plot_abundance(varargin)
function plot_abundance(f,c,outname,tit);

%    varargin: list of input files to plot on abundance plot

%c={'bo','ro','go','co','ko','ks','bs','rs','gs','cs'};
%c={'bo','ro','go','co','ko'};
%c={'ks','bs','rs','gs','cs'};


sample_name={'TCR6';'TCR7';'TCR22';'TCR23';'TCR28';'TCR29';'TCR12';'TCR13';'TCR14';'TCR15';'TCR8';'TCR9'};
sample_type={'LGG1';'LGG1';'LGG2';'LGG2';'LGG4';'LGG4';'GBM2';'GBM2';'GMB1';'GBM1';'normal';'normal'};

%{
sample_name={'TCR40'
    'TCR41'
    'TCR42'
    'TCR43'
    'TCR6'
    'TCR7'
    'SIMS0101-PBMC.tsv.VJ.productive'
    'SIMS0101-CD4+CD25+.tsv.VJ.productive'
    'SIMS0101-CD4-CD8.tsv.VJ.productive'
    'SIMS0101-CD8+CD56+.tsv.VJ.productive'
    'SIMS0101-TUMOR.tsv.VJ.productive'};
sample_type={'NKT (iRepertoire)','CD8 (iRepertoire)','Treg (iRepertoire)','CD4 (iRepertoire)','LGG brain (iRepertoire)','LGG blood (iRepertoire)','LGG blood (Immunoseq)','Treg (Immunoseq)','NKT (Immunoseq)','CD8 (Immunoseq)','LGG brain (Immunoseq)'};
%}
leg={};
fig=figure;
k=1;
for i=1:length(f)
    [pathstr,name,ext] = fileparts(f{i});
    fin=fopen(f{i},'r');
    h=fgetl(fin);
    data=textscan(fin,'%s%d%s%s%s%d');
    cdr3=data{5};
    counts=data{2};
    fclose(fin);
    
    cc=find(counts>1);
    cdr3=cdr3(cc);
    counts=counts(cc);
    

    [cdr3,ind]=sort(cdr3); % sort sequences
    counts=counts(ind); % sort counts
    [u_aa,uA,uS]=unique(cdr3); % unique amino acid sequences and indices

    % summed counts
    aacounts=zeros(length(u_aa),1);
    for j=1:length(uS)
	aacounts(uS(j))=aacounts(uS(j))+counts(j);
    end

    % count number of distinct CDR3s with that clone size
    uc=unique(aacounts);
    nunq=histc(aacounts,uc);
    name_brev=strtok(name,'_');
%	name_brev
    annot=strcmp(name_brev,sample_name);
    leg{i}=sample_type{annot}
 %   if k<=n1
%	c=c1;
    %end
    %if k==(n1+1);
    %  c=c2;
    %   k=1;
    %end 

    loglog(uc,nunq,c{i},'MarkerSize',6)
    k=k+1;
    hold on
end
    title(tit,'FontSize',16)
    xlabel('copy #','FontSize',16)
    ylabel('unique CDR3s','FontSize',16)
    set(gca,'FontSize',14) 
    legend(leg,'FontSize',14)

    outdir='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/analysis/plots/abundance/';
    outpath=[outdir,outname];
    print(fig,'-dpdf','-r300',[outpath,'.pdf']);

