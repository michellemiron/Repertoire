function get_entropy(s,datafile,epath,vjpath,varargin)

%{
cases
0 - cdr3
1 - vj both
2 - v only
3 - j only
%}

% read in file
fin=fopen(datafile);
header = fgetl(fin);
data=textscan(fin,'%s%f%s%s%s%f');

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

% split datafile into parts
[pathstr,name,ext] = fileparts(datafile);
% process arguments

if nargin==4
	cutoff=length(unique(cdr3));
	label='all';
elseif nargin==5
	cutoff=varargin{1};
	label=num2str(cutoff);
else
	disp('error: incorrect number of arguments')
	exit;
end


% merge identical elements and compute entropy
switch s
	case 0 % cdr3
	   Htype='Hcdr3';
	   % sort CDR3s alphabetically
	   [cdr3,ind]=sort(cdr3); % sort sequences
	   counts=counts(ind); % sort counts
	   [u_aa,uA,uS]=unique(cdr3); % unique amino acid sequences and indices
	   totalnum=length(u_aa(1:cutoff));
            % summed counts
            aacounts=zeros(length(u_aa),1);        
            for j=1:length(uS)
                aacounts(uS(j))=aacounts(uS(j))+counts(j);
            end
		% sort by increading counts
		[aacounts,k]=sort(aacounts,'descend');
		u_aa=u_aa(k);

		% compute entropy
		f=aacounts(1:cutoff)/sum(aacounts(1:cutoff));
		H=entropy(f);
	
	case 1 % vj
	    
	    Htype='Hvj';
		keyboard
	    % get unique cdr3s in counts order 
	    [u_aa,uA,uS]=unique(cdr3,'first');
	    u_aa_sorted=cdr3(sort(uA));
	
	    % determine cutoff cdr3 and select corresponding VJ with counts
	    cutoffcdr3=u_aa_sorted(cutoff);
	    i=find(strcmp(cutoffcdr3,cdr3));
	    vj=vj(1:i(end));
	    counts=counts(1:i(end));
	    
	    [u_vj,uA,uS]=unique(vj);
	   
	    totalnum=length(u_vj);
	    vjcounts=zeros(length(u_vj),1);
	    for j=1:length(uS)
                vjcounts(uS(j))=vjcounts(uS(j))+counts(j);
            end
	    
	    % get all V segments
	    Vin=fopen([pathstr,'/allV'],'r'); 
	    Vsegs=textscan(Vin,'%s'); Vsegs=Vsegs{1};
	    % get all J segments
	    Jin=fopen([pathstr,'/allJ'],'r');
	    Jsegs=textscan(Jin,'%s'); Jsegs=Jsegs{1};

	    fclose(Vin); fclose(Jin);
	   
	    % convert counts to frequencies
	    f=vjcounts/sum(vjcounts);
	    
	   % save VJ frequencies in a table
	    vjpath=[vjpath,'/',name,'.',label,'.VJuse.txt'];
	    fout=fopen(vjpath,'w');
	    for i=1:length(Vsegs)
		vi=Vsegs{i};
		for j=1:length(Jsegs)
		    ji=Jsegs{j};
		    ind=find(strcmp([vi,'.',ji],u_vj));
		    if(~isempty(ind))
			fprintf(fout,[vi,'\t',ji,'\t',num2str(f(ind)),'\n']); % write in counts where VJ pairing is observed
		    else
			fprintf(fout,[vi,'\t',ji,'\t',num2str(0),'\n']); % write in 0 where VJ pairing is not observed
		    end
		end 
	    end
	   fclose(fout);
	   H=entropy(f); % compute entropy
	    

	case 2 % v
	    
	    Htype='Hv';
	    [u_aa,uA,uS]=unique(cdr3,'first');
	    u_aa_sorted=cdr3(sort(uA));
	
	    % determine cutoff cdr3 and select corresponding VJ with counts
	    cutoffcdr3=u_aa_sorted(cutoff);
	    i=find(strcmp(cutoffcdr3,cdr3));
	    vseg=vseg(1:i(end));
	    counts=counts(1:i(end));
	    
	    [u_v,uA,uS]=unique(vseg);
	    totalnum=length(u_v);
	    vcounts=zeros(length(u_v),1);
	    for j=1:length(uS)
                vcounts(uS(j))=vcounts(uS(j))+counts(j);
            end

            % get all V segments
            Vin=fopen([pathstr,'/allV'],'r');
            Vsegs=textscan(Vin,'%s'); Vsegs=Vsegs{1};

	    % convert counts to frequencies
	    f=vcounts/sum(vcounts);

	    vpath=[vjpath,'/',name,'.',label,'.Vuse.txt'];
	    fout=fopen(vpath,'w');
            for i=1:length(Vsegs) 
		vi=Vsegs{i};
                ind=find(strcmp(vi,u_v));
                    if(~isempty(ind))
                        fprintf(fout,[vi,'\t',num2str(f(ind)),'\n']); % write in counts where VJ pairing is observed
                    else
                        fprintf(fout,[vi,'\t',num2str(0),'\n']); % write in 0 where VJ pairing is not observed
                    end
            end
            fclose(fout);
	    H=entropy(f); % compute entropy

	case 3 % j
	    Htype='Hj';
	    [u_aa,uA,uS]=unique(cdr3,'first');
	    u_aa_sorted=cdr3(sort(uA));
	
	    % determine cutoff cdr3 and select corresponding VJ with counts
	    cutoffcdr3=u_aa_sorted(cutoff);
	    i=find(strcmp(cutoffcdr3,cdr3));
	    jseg=jseg(1:i(end));
	    counts=counts(1:i(end));
	     
	    [u_j,uA,uS]=unique(jseg);
	    totalnum=length(u_j);
	    jcounts=zeros(length(u_j),1);
	    for j=1:length(uS)
                jcounts(uS(j))=jcounts(uS(j))+counts(j);
            end

            % get all J segments
            Jin=fopen([pathstr,'/allJ'],'r');
            Jsegs=textscan(Jin,'%s'); Jsegs=Jsegs{1};
	   
  	    % convert counts to frequencies
  	    f=jcounts/sum(jcounts);

            jpath=[vjpath,'/',name,'.',label,'.Juse.txt'];
            fout=fopen(jpath,'w');
            for i=1:length(Jsegs)
                ji=Jsegs{i};
                ind=find(strcmp(ji,u_j));
                    if(~isempty(ind))
                        fprintf(fout,[ji,'\t',num2str(f(ind)),'\n']); % write in counts where VJ pairing is observed
                    else
                        fprintf(fout,[ji,'\t',num2str(0),'\n']); % write in 0 where VJ pairing is not observed
                    end
            end
	    fclose(fout);
	    H=entropy(f);

	otherwise
		disp('error: unspecified case')
		exit;
end

outpath=[epath,'/',name,'.',label,'.',Htype];
fout=fopen(outpath,'w');
fprintf(fout,[num2str(totalnum),'\t', num2str(H),'\n']);
fclose(fout);
