%{ 
% iRep vs Immunoseq
f1='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR40_both.trimed.ready.bam.Bchain.productive.tsv';
f2='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR41_both.trimed.ready.bam.Bchain.productive.tsv';
f3='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR42_both.trimed.ready.bam.Bchain.productive.tsv';
f4='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR43_both.trimed.ready.bam.Bchain.productive.tsv';
f5='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR6_S1_L001_R1_001.Bchain.productive.tsv';
f6='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR7_S1_L001_R1_001.Bchain.productive.tsv';
f7='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/immunoseq/SIMS0101-PBMC.tsv.VJ.productive.tsv';
f8='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/immunoseq/SIMS0101-CD4+CD25+.tsv.VJ.productive.tsv';
f9='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/immunoseq/SIMS0101-CD4-CD8.tsv.VJ.productive.tsv';
f10='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/immunoseq/SIMS0101-CD8+CD56+.tsv.VJ.productive.tsv';
f11='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/immunoseq/SIMS0101-TUMOR.tsv.VJ.productive.tsv';
f={f1;f2;f3;f4;f5;f6;f9;f10;f8;f11;f7}; %   f7;f8;f9;f10;f11};
%c={'bo';'ro';'go';'mo';'co';'ko';'ks';'gs';'bs';'rs';'cs'};
c={'bo';'ro';'go';'mo';'co';'ko';'bs';'rs';'gs';'cs';'ks'};

inds=7:11;
f=f(inds);
c=c(inds);
plot_abundance(f,c,'TCR40-43_abundance_Immunoseq');
%}

%{
% A chain
f1='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRA/Jenny_out/TCR6_S1_L001_R1_001.Achain.productive.tsv';
g1='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRA/Jenny_out/TCR7_S1_L001_R1_001.Achain.productive.tsv';
f2='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRA/Jenny_out/TCR22_both.Achain.productive.tsv';
g2='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRA/Jenny_out/TCR23_both.Achain.productive.tsv';
f3='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRA/Jenny_out/TCR28_both.Achain.productive.tsv';
g3='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRA/Jenny_out/TCR29_both.Achain.productive.tsv';
f4='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRA/Jenny_out/TCR12_S1_L001_R1_001.Achain.productive.tsv';
g4='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRA/Jenny_out/TCR13_S1_L001_R1_001.Achain.productive.tsv';
f5='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRA/Jenny_out/TCR14_S1_L001_R1_001.Achain.productive.tsv';
g5='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRA/Jenny_out/TCR15_S1_L001_R1_001.Achain.productive.tsv';
f6='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRA/Jenny_out/TCR8_S1_L001_R1_001.Achain.productive.tsv';

bl={f1;f2;f3;f4;f5;f6};
br={g1;g2;g3;g4;g5};
c_bl={'bo';'ro';'go';'mo';'co';'ko'};
c_br={'bs';'rs';'gs';'ms';'cs'};

plot_abundance(br,c_br,'Achain_abundance_brain','Human Brain - A chain');
%}

% B chain
g1='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR6_S1_L001_R1_001.Bchain.productive.tsv';
f1='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR7_S1_L001_R1_001.Bchain.productive.tsv';
g2='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR22_both.Bchain.productive.tsv';
f2='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR23_both.Bchain.productive.tsv';
g3='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR28_both.Bchain.productive.tsv';
f3='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR29_both.Bchain.productive.tsv';
g4='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR12_S1_L001_R1_001.Bchain.productive.tsv';
f4='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR13_S1_L001_R1_001.Bchain.productive.tsv';
g5='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR14_S1_L001_R1_001.Bchain.productive.tsv';
f5='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR15_S1_L001_R1_001.Bchain.productive.tsv';
f6='/ifs/scratch/c2b2/ys_lab/bg2178/projects/tcr/Sims/dat/processedbams/human/TCRB/Jenny_out/TCR9_.Bchain.productive.tsv';

bl={f1;f2;f3;f4;f5;f6};
br={g1;g2;g3;g4;g5};
c_bl={'bo';'ro';'go';'mo';'co';'ko'};
c_br={'bs';'rs';'gs';'ms';'cs'};

plot_abundance(br,c_br,'Bchain_abundance_brain','Human Brain - B chain');
