function get_entropy_byVJ(datafile)

% read in file
fin=fopen(datafile);
header=fgetl(fin);
data=textscan(fin,'%s%f%s%s%s%f');
[pathstr,name,ext] = fileparts(datafile);

vj=data{1};
counts=data{2};
cdr3=data{5};
vseg=data{3};
jseg=data{4};
cdr3len=data{6};
fclose(fin);

% remove clones of size 1
cc=find(counts>1);
vj=vj(cc);
cdr3=cdr3(cc);
vseg=vseg(cc);
jseg=jseg(cc);
cdr3len=cdr3len(cc);
counts=counts(cc);

% Full CDR3 entropy
[cdr3,ind]=sort(cdr3); % sort sequences
counts=counts(ind); % sort counts
[u_aa,uA,uS]=unique(cdr3); % unique amino acid sequences and indices

% summed counts
aacounts=zeros(length(u_aa),1);
for j=1:length(uS)
    aacounts(uS(j))=aacounts(uS(j))+counts(j);
end

% compute entropy
f=aacounts/sum(aacounts);
Hcdr3=entropy(f);


% Individual VJ entropies
uvj=unique(vj); % unique VJs
L=length(uvj);
H=zeros(L,1);
num_vj=zeros(L,1);
for i=1:L
	s=strcmp(uvj(i),vj); % CDR3s corresponding to a given VJ
	num_vj(i)=sum(s); % number associated CDR3a
	aa=cdr3(s); % amino acid sequences
	vals=counts(s); % values
	f=vals/sum(vals);
	maxH=-log2(1/length(vals));
	H(i)=entropy(f);%/maxH; % compute entropy
	H(isnan(H))=0;
end

% plot

outdir='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/analysis/plots/entropies';

% Change sorting pattern here %
%%%%%%%%%%%%

%[H,inds]=sort(H,'descend');
%num_vj=num_vj(inds);
%outdir=[outdir,'/',name,'_','entropy_byVJ_Hsort'];

[num_vj,inds]=sort(num_vj,'descend');
H=H(inds);
outdir=[outdir,'/',name,'_','entropy_byVJ_freqsort'];
%%%%%%%%%%%%


uvj=uvj(inds);

fig=figure;

set(gcf,'PaperUnits', 'inches','PaperPosition',[0.1,0.1,9,11]);

%plot(1:L,H./log2(num_vj),'o')
plot(1:L,log2(num_vj),'ro');%,'Linewidth',2)
hold on
bar(1:L,H);


%bar(1:L,H./log2(num_vj))

title('Entropy by VJ');

% adjust x-axis
%{
Xl = [1:L];
set(gca,'XTick',Xl);%,'XLim',[1, length(fname)]);
set(gca,'XTickLabel','')
set(gca,'FontSize',18)

ax = axis; % current axis limits
axis(axis); % set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4); % y-axis limits

% place the text labels
t = text(Xl,Yl(1)*ones(1,length(Xl)),uvj,'FontSize',14,'FontWeight','bold');
set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
'Rotation',45);
%}

xlim([1,L])
legend('max entropy (log scaled frequencies)', 'VJ entropy')
xlabel('VJ');
ylabel('H');


grid on
box on

print(fig,'-dpdf','-r300',[outdir,'.pdf'])
