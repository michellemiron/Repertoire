function plot_segs(fpath)
%
% Inputs:
%	fpath: directory path to file
%	cutoff: switch. 0 for no cutoff, 1 for cutoff currently preset to 1000
%

addpath '/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/analysis/mplot' % for colors

% load file and collect data
fin=fopen(fpath,'r');
header=fgetl(fin);
header=regexp(deblank(header),'\t','split');
hl=length(header)-3;

data=textscan(fin,['%s%s%s', repmat('%f',1,hl)]','delimiter', '\t', 'EmptyValue', -Inf);
fname=data{1};
tissue=data{2};
sample=data{3};
vals=cell2mat(data(4:end));

[pathstr,name,ext] = fileparts(fpath);
[name,seg]=strtok(name,'.');
[cutoff,seg]=strtok(seg,'.');
seg=seg(2:end);
if strcmp('humanAchain',name)
	tit=[seg, ' cassette for human TCR-alpha' ];
elseif strcmp('humanBchain',name)
	tit=[seg,' cassette for human TCR-beta'];
elseif strcmp('mouseAchain',name)
	tit=[seg ' cassette for mouse TCR-alpha'];
elseif strcmp('mouseBchain',name)
	tit=[seg, ' cassette for mouse TCR-beta'];
else
	'what is this file?'
%	exit()
end

tit=[seg, ' cassette for iRepertoire and Immunoseq'];

% separate brain and blood data
inds_bl=find(strcmp('blood',tissue));
inds_br=find(strcmp('brain',tissue));
fname_bl=fname(inds_bl); samples_bl=sample(inds_bl); vals_bl=vals(inds_bl,:);
fname_br=fname(inds_br); samples_br=sample(inds_br); vals_br=vals(inds_br,:);

% pair brain-blood samples and construct cell containing blood in first entry and brain in second
% where there is no pairing for brain, put in 0s as filler.
groupedcell={};
for i=1:length(fname_bl)
s=samples_bl(i);
if sum(strcmp(s,samples_br))
	groupedcell{1}(i,:)=vals_bl(i,:);
	groupedcell{2}(i,:)=vals_br(strcmp(s,samples_br),:);
else
	groupedcell{1}(i,:)=vals_bl(i,:);
end

end

% plot
fig=figure;
%set(gcf,'PaperUnits', 'inches','PaperPosition',[0.25,0.25,9,9],'PaperSize',[10,10]);
set(gcf,'PaperUnits', 'inches','PaperPosition',[0.25,0.25,14,14],'PaperSize',[15,15]); % for mouse V
bar(groupedcell{1},'stacked','BarWidth',0.3)
hold on
Xl = [1:length(samples_bl)];
bar(groupedcell{2},'stacked','BarWidth',0.3,'Xdata', Xl+0.3);
xlim([0.5,Xl(2)+5])
ylim([0,1.05])
title(tit,'fontsize',14);
ylabel('frequency','FontSize',12);
colormap(distinguishable_colors(length(header)-3));

gridLegend(header(4:end),2,'location','BestOutside','Fontsize',10);
%gridLegend(header(4:end),2,'location','BestOutside','Fontsize',5); % for mouse V

% adjust x-axis
set(gca,'XTick',Xl);%,'XLim',[1, length(fname)]);
set(gca,'XTickLabel','')
set(gca,'FontSize',10)

ax = axis; % current axis limits
axis(axis); % set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4); % y-axis limits

% place the text labels
t = text(Xl,Yl(1)*ones(1,length(Xl)),samples_bl,'FontSize',10,'FontWeight','bold');
set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
'Rotation',45);
u=text(Xl,Yl(2)*ones(1,length(Xl))-0.05,'iRepertoire','FontSize',10);
set(u,'Rotation',45','color','r');
Xll=[Xl(1:3), Xl(5:end)];
%Xll=Xl(2:end);

v=text(Xll+0.3,Yl(2)*ones(1,length(Xll))-0.05,'Immunoseq','FontSize',10);
set(v,'Rotation',45','color','b');


box on

outdir='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/src2/analysis/plots/vj';
outdir=[outdir,'/',name,'.',cutoff,'.',seg];
print(fig,'-dpdf','-r300',[outdir,'.pdf'])
