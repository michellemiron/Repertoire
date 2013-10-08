function get_entropy3(s,datafile,varargin)

% Same as get_entropy2, but doesn't write to file and takes fewer inputs for eaier commandline run


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

if nargin==2
	cutoff=length(unique(cdr3));
	label='all';
elseif nargin==3
	cutoff=varargin{1};
	label=num2str(cutoff);
else
	disp('error: incorrect number of arguments')
end


% merge identical elements and compute entropy
switch s
	case 0 % cdr3
	   Htype='Hcdr3';
	   
	   % sort alphabetically 
	   [cdr3,ind]=sort(cdr3); % sort sequences
	   counts=counts(ind); % sort counts
	
	   [u_aa,uA,uS]=unique(cdr3); % unique amino acid sequences and indices
	   totalnum=length(u_aa);
            % summed counts
            aacounts=zeros(length(u_aa),1);        
            for j=1:length(uS)
                aacounts(uS(j))=aacounts(uS(j))+counts(j);
            end
		% sort by increasing counts
		[aacounts,k]=sort(aacounts,'descend');
		u_aa=u_aa(k);
		
		% apply cutoff
                u_aa=u_aa(1:cutoff);
                aacounts=aacounts(1:cutoff);
		
		% compute entropy
		f=aacounts/sum(aacounts);
		H=entropy(f)
	
	case 1 % vj
	    
	    Htype='Hvj';
	    % get unique cdr3s in counts order 
	    %[u_aa,uA,uS]=unique(cdr3,'first');
	    %u_aa_sorted=cdr3(sort(uA));
	
	    % determine cutoff cdr3 and select corresponding VJ with counts
	    %cutoffcdr3=u_aa_sorted(cutoff);
	    %i=find(strcmp(cutoffcdr3,cdr3));
	    [u_vj,uA,uS]=unique(vj);
	     totalnum=length(u_vj);
	    vjcounts=zeros(length(u_vj),1);
	    for j=1:length(uS)
                vjcounts(uS(j))=vjcounts(uS(j))+counts(j);
            end
	   keyboard
	    vj=vj(1:cutoff);
            counts=counts(1:cutoff);
	    % convert counts to frequencies
	    f=vjcounts/sum(vjcounts);
	    
	   H=entropy(f) % compute entropy
	    

	case 2 % v
	    
	    Htype='Hv';
	    %[u_aa,uA,uS]=unique(cdr3,'first');
	    %u_aa_sorted=cdr3(sort(uA));
	
	    % determine cutoff cdr3 and select corresponding VJ with counts
	    %cutoffcdr3=u_aa_sorted(cutoff);
	    %i=find(strcmp(cutoffcdr3,cdr3));
	    vseg=vseg(1:cutoff);
	    counts=counts(1:cutoff);
	    
	    [u_v,uA,uS]=unique(vseg);
	    totalnum=length(u_v);
	    vcounts=zeros(length(u_v),1);
	    for j=1:length(uS)
                vcounts(uS(j))=vcounts(uS(j))+counts(j);
            end

	    % convert counts to frequencies
	    f=vcounts/sum(vcounts);

	    H=entropy(f) % compute entropy

	case 3 % j
	    Htype='Hj';
	    %[u_aa,uA,uS]=unique(cdr3,'first');
	    %u_aa_sorted=cdr3(sort(uA));
	
	    % determine cutoff cdr3 and select corresponding VJ with counts
	    %cutoffcdr3=u_aa_sorted(cutoff);
	    %i=find(strcmp(cutoffcdr3,cdr3));
	    jseg=jseg(1:cutoff);
	    counts=counts(1:cutoff);
	     
	    [u_j,uA,uS]=unique(jseg);
	    totalnum=length(u_j);
	    jcounts=zeros(length(u_j),1);
	    for j=1:length(uS)
                jcounts(uS(j))=jcounts(uS(j))+counts(j);
            end

  	    f=jcounts/sum(jcounts);

	    H=entropy(f)

	otherwise
		disp('error: unspecified case')
		exit;
end

