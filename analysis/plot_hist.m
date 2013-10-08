function plot_hist(f_bl,f_br,type,chain,order)

fin1=fopen(f_bl);
fin2=fopen(f_br);

% headers
h=fgetl(fin1);
fgetl(fin2);

% load files
data_bl=textscan(fin1,'%s%d%s%s%s%d');
data_br=textscan(fin2,'%s%d%s%s%s%d');

fclose(fin1); fclose(fin2);

% extract data
cdr3_bl=data_bl{5};
vals_bl=data_bl{2};
cdr3_br=data_br{5};
vals_br=data_br{2};

% remove clones of size 1
cc=find(vals_bl>1);
cdr3_bl=cdr3_bl(cc);
vals_bl=vals_bl(cc);

cc=find(vals_br>1);
cdr3_br=cdr3_br(cc);
vals_br=vals_br(cc);


% add counts for identical cdr3s

[cdr3_bl,ind]=sort(cdr3_bl); % sort sequences
vals_bl=vals_bl(ind); % sort counts
[u_bl,uA,uS]=unique(cdr3_bl); % unique amino acid sequences and indices

% summed counts
blcounts=zeros(length(u_bl),1);
for j=1:length(uS)
    blcounts(uS(j))=blcounts(uS(j))+vals_bl(j);
end

[cdr3_br,ind]=sort(cdr3_br); % sort sequences
vals_br=vals_br(ind); % sort counts
[u_br,uA,uS]=unique(cdr3_br); % unique amino acid sequences and indices

% summed counts
brcounts=zeros(length(u_br),1);
for j=1:length(uS)
    brcounts(uS(j))=brcounts(uS(j))+vals_br(j);
end

% include any CDR3s that are in one sample and not in another (give 0 count)
blonly=setdiff(u_bl,u_br); % only in blood sample
bronly=setdiff(u_br,u_bl); % only in brain sample

u_bl=[u_bl; bronly]; blcounts=[blcounts; zeros(length(bronly),1)];
u_br=[u_br; blonly]; brcounts=[brcounts; zeros(length(blonly),1)];

% sort cdr3s alphabetically (u_bl and u_br should be in identical order)
[u_bl,ind]=sort(u_bl);
blcounts=blcounts(ind);
[u_br,ind]=sort(u_br);
brcounts=brcounts(ind);

% sort by specified order
switch order
case 0    % sort all according to blood counts
[blcounts,ind]=sort(blcounts,'descend');
u_bl=u_bl(ind);
u_br=u_br(ind);
brcounts=brcounts(ind);
blfreq=blcounts/sum(blcounts);
brfreq=brcounts/sum(brcounts);

% write to file
outdir='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/analysis/counts';
outfile=[outdir,'/',type,'_bloodorder_',chain,'chain.tsv'];
fout=fopen(outfile,'w');
fprintf(fout,'blood_CDR3\tblood_counts\tblood_frequency\tbrain_CDR3\tbrain_counts\tbrain_frequency\n');
for i=1:length(u_bl)
		fprintf(fout,[u_bl{i},'\t',num2str(blcounts(i)),'\t',num2str(blfreq(i)),'\t',u_br{i},'\t',num2str(brcounts(i)),'\t',num2str(brfreq(i)),'\n']);
end

case 1
[brcounts,ind]=sort(brcounts,'descend'); % sort all according to brain counts
u_bl=u_bl(ind);
u_br=u_br(ind);
blcounts=blcounts(ind);
blfreq=blcounts/sum(blcounts);
brfreq=brcounts/sum(brcounts);

% write to file
outdir='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/analysis/counts';
outfile=[outdir,'/',type,'_brainorder_',chain,'chain.tsv'];
fout=fopen(outfile,'w');
fprintf(fout,'blood_CDR3\tblood_counts\tblood_frequency\tbrain_CDR3\tbrain_counts\tbrain_frequency\n');
        for i=1:length(u_bl)
                fprintf(fout,[u_bl{i},'\t',num2str(blcounts(i)),'\t',num2str(blfreq(i)),'\t',u_br{i},'\t',num2str(brcounts(i)),'\t',num2str(brfreq(i)),'\n']);
	end
end
fclose(fout);


%keyboard
% plot
fig=figure;

switch order
case 0
	bar(brfreq,100,'r','Edgecolor','None'); hold on; bar(blfreq,70,'b','Edgecolor','None');
	legend({'brain','blood'},'FontSize',16);
	title(['Blood-brain histogram for ', type,', ',chain,' chain (blood order)'] ,'fontsize',20);
case 1 
        bar(blfreq,100,'b','Edgecolor','None'); hold on; bar(brfreq,70,'r','Edgecolor','None');
	legend({'blood','brain'},'FontSize',16);
        title(['Blood-brain histogram for ', type,', ',chain,' chain (brain order)'] ,'fontsize',20);
	bar(blfreq,70,'b','Edgecolor','None'); hold on; bar(brfreq,70,'r','Edgecolor','None');
end

set(gcf,'PaperUnits', 'inches','PaperPosition',[0.1,0.1,8,8]);

% xlabel('samples')

ylabel('frequency','FontSize',14);

% adjust x-axis
xlabel('ordered CDR3s','FontSize',14)

ax = axis; % current axis limits
axis(axis); % set the axis limit modes (e.g. XLimMode) to manual
set(gca,'XTick',[])

Yl = ax(3:4); % y-axis limits

xlim([0,length(u_bl)])
%xlim([0,20000])
ylim([0 0.07]);
%grid on
box on

outdir='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/analysis/plots/histograms';
switch order
case 0
	outfile=[outdir,'/',type,'_bloodorder_',chain,'chain'];
case 1
        outfile=[outdir,'/',type,'_brainorder_',chain,'chain'];
end

set(gcf,'Renderer','painters')
print(fig,'-dpdf','-r300',[outfile,'.pdf'])

