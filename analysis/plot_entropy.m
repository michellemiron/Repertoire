function plot_entropy(fpath,outdir,bb,cutoff)
% Inputs
%     fpath: directory path to input file
%     outdir: output directory
%     bb: tissue type (brain or blood)
%     cutoff: switch. 0 for no cutoff, 1 for cutoff currently preset to 1000
%

% read file
fin=fopen(fpath,'r');
header=fgetl(fin);
data=textscan(fin,'%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter', '\t', 'EmptyValue', -Inf);

fname=data{1};
tissue=data{2};
sample=data{3};
ncdr3=data{4};
Hcdr3=data{5};
Hcdr3_1000=data{6};
nvj=data{7};
Hvj=data{8};
nvj_1000=data{9};
Hvj_1000=data{10};
nv=data{11};
Hv=data{12};
nv_1000=data{13};
Hv_1000=data{14};
nj=data{15};
Hj=data{16};
nj_1000=data{17};
Hj_1000=data{18};

fclose(fin);
[pathstr,name,ext] = fileparts(fpath);
name=strtok(name,'.');

if strcmp('humanAchain',name)
	tit=['TCR-alpha for human ', bb, ' samples' ];
elseif strcmp('humanBchain',name)
	tit=['TCR-beta for human ', bb, ' samples' ];
elseif strcmp('mouseAchain',name)
	tit=['TCR-alpha for mouse ', bb, ' samples' ];
elseif strcmp('mouseBchain',name)
	tit=['TCR-beta for mouse ', bb, ' samples' ];
else
	'what is this file?'
%	exit()
end

% no cutoff or yes cutoff (set to 1000 for now)
switch cutoff
	case 0
		maxentropy=-log2(1./ncdr3);
		a='all';
	case 1
		maxentropy=repmat(-log2(1/1000),length(fname),1);
		Hcdr3=Hcdr3_1000;
		nvj=nvj_1000;
		Hvj=Hvj_1000;
		nv=nv_1000;
		Hv=Hv_1000;
		nj=nj_1000;
		Hj=Hj_1000;
		a='1000';
end

% plot
fig=figure;
set(gcf,'PaperUnits', 'inches','PaperPosition',[0.1,0.1,9,11]);

inds=find(strcmp(bb,tissue));

H1=Hcdr3(inds);
H2=Hvj(inds)./H1;
H3=Hv(inds)./H1;
H4=Hj(inds)./H1;
H1=H1./H1;
%H1=Hcdr3(inds)./maxentropy(inds); H1(isnan(H1))=0; H1(isinf(H1))=0;
%H2=Hvj(inds)./log2(nvj(inds)); H2(isnan(H2))=0; H2(isinf(H2))=0;
%H3=Hv(inds)./log2(nv(inds)); H3(isnan(H3))=0; H3(isinf(H3))=0;
%H4=Hj(inds)./log2(nj(inds)); H4(isnan(H4))=0; H4(isinf(H4))=0;

bar([H1,H2,H3,H4],'grouped');
%hold on
%plot(1:length(ncdr3(inds)),maxentropy(inds),'*');
legend({'Hcdr3';'Hvj';'Hv';'Hj'},'location','BestOutside','FontSize',20);
title(tit,'fontsize',22);
% xlabel('samples')

ylabel('normalized entropy','FontSize',20);

% adjust x-axis
Xl = [1:length(sample(inds))];
set(gca,'XTick',Xl);%,'XLim',[1, length(fname)]);
set(gca,'XTickLabel','')
set(gca,'FontSize',18)

ax = axis; % current axis limits
axis(axis); % set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4); % y-axis limits

% place the text labels
t = text(Xl,Yl(1)*ones(1,length(Xl)),sample(inds),'FontSize',14,'FontWeight','bold');
set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
'Rotation',45);

ylim([0 1.05]);
grid on
box on

outdir=[outdir,'/',name,'_','entropies','_',a,'_',bb];
print(fig,'-dpdf','-r300',[outdir,'.pdf'])

%saveas(fig,outdir,'png')
