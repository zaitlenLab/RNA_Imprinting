% This code estimates the parameters of the three beta distributions which
% represent the balanced, imbalanced and imprinted expression classes

% !! replace with the directory into which imprinting_code was extracted
refdir = '/Users/kthefrog/path_to_dir/imprinting_code';

%=====================
% data files and flags
%=====================

% directories of code files
addpath([refdir '/parameter_inference']);
addpath([refdir '/common']);

% directory of data files
dprefix = [refdir '/data'];
auxprefix = [refdir '/aux_files'];

% dataset used for inference
datafile = 'GD462.ASE.COV16.ANNOTPLUS.SINFO.txt.LCL';

% (individual,gene) pairs for which NMD event is predicted
% these will be omitted from the analysis
nmdfile = [auxprefix '/' 'GD_ind_gene_NMD.txt'];
% a list of types (as opposed to imputed) SNPs
typedfile = [auxprefix '/' 'GEUVADIS.PH1PH2_465.GTd.snps.txt'];
% a list of individual from whom phasing in not given
gdphasefile = [auxprefix '/' 'imputed_samples_phase_error.txt'];

%==========
% READ DATA
%==========

disp(['reading files for ' datafile]);

%========= load count data (reference, total counts) for informative sites

load([dprefix '/' datafile '.mat'],'ref','n');
snpn = size(n,2);

%========= load gene annotations

load([dprefix '/' datafile '.mat'],'sgen');
[u_gen i2o_gen i2u_gen] = unique(sgen);
genn = length(u_gen);

%========= load subject ids

id = textread([dprefix '/' datafile '.sub'],'%s');
[u_id i2o_id i2u_id] = unique(id);
indn = length(u_id);

%========= process NMD data

[nmd_id nmd_gen] = textread(nmdfile,'%s%s');
nmd = false(size(n));
for i = (1:length(nmd_id)),
    thisid = find(strcmp(nmd_id(i),u_id),1,'first');
    thisgen = find(strcmp(nmd_gen(i),u_gen),1,'first');
    if (~isempty(thisgen)),
        nmd(thisid,i2u_gen==thisgen) = true;

    end
end

%========= FIGURE OUT WHICH SNPS ARE IMPUTED

rsfile = [dprefix '/' datafile '.rs'];
rs = textread(rsfile,'%s');
[u_rs i2o_rs i2u_rs] = unique(rs);
posfile = [dprefix '/' datafile '.pos'];
pos = load(posfile);
chrfile = [dprefix '/' datafile '.chr'];
chr = load(chrfile);

[chrpos ind] = sortrows([chr(i2o_rs) pos(i2o_rs)],[1 2]);
u_rs = u_rs(ind);
typed = textread(typedfile,'%s');
imputed = ~ismember(u_rs,typed);
clear pos chr ind chrpos; 

%===========================
% FILTER (minimal filtering)
%===========================

disp('filtering SNPs');

remove = false(snpn,1);

%========= HWE FILTERING

load([dprefix '/' datafile '.hwe.mat'],'hwep');
hweremove = (0<=hwep & hwep<1e-2)';
remove = (remove|hweremove);

remove = (remove|i2u_gen==1);

% clean removed SNPs
ref = ref(:,~remove);
n = n(:,~remove);
i2u_gen = i2u_gen(~remove);
imputed = imputed(~remove);
u_rs = u_rs(~remove);

% clean NMD sites
nmd = nmd(:,~remove);
n(nmd) = 0;
ref(nmd) = 0;
clear nmd;

%========= switch from matrix to vector for efficiency

indn = size(ref,1);
snpn = size(ref,2);
nonan = (n>0 & n<300);
ref = ref(nonan); % ref counts of all relevant sites
alt = n(nonan)-ref; % alt counts of all relevant sites
clear n;
tmp = repmat((i2u_gen)',[indn 1]);
gene = tmp(nonan);
tmp = repmat((1:indn)',[1 snpn]);
ind = tmp(nonan);
tmp = repmat(1:snpn,[indn 1]);
snp = tmp(nonan);

%========= obtain initial values for beta params based

% per gene, compute the mean asbolute deviation from 0.5

genef = zeros(genn,1);
for i = (1:genn),
   genef(i) = mean(abs(ref(gene==i)./(ref(gene==i)+alt(gene==i))-0.5)); 
end
q = quantile(genef,0:0.1:1);

% use quantiles of genef to initialize beta param estimates

params = zeros(3,2);
starts = [2 2]; % initial values passed to optimization function
smallbnd = 1e-5;
options = optimset('algorithm','interior-point');

% estimate initial parameters for balanced class using 1st decile
cursites = (q(1)<=genef & genef<=q(2));
cursites = cursites(gene);
fun = @(x)betamle_1(x,ref(cursites),alt(cursites));
params(1,:) = fmincon(fun,starts,[],[],[],[], ...
        [smallbnd smallbnd],[Inf Inf],[],options);

% estimate initial parameters for imbalanced class using 9th decile
cursites = (q(9)<=genef & genef<=q(10));
cursites = cursites(gene);
fun = @(x)betamle_23(x,ref(cursites),alt(cursites));
params(2,:)  = fmincon(fun,starts,[],[],[],[], ...
        [smallbnd smallbnd],[Inf Inf],[],options);

% estimate initial parameters for imprinted class using 10th decile
cursites = (q(10)<=genef & genef<=q(11));
cursites = cursites(gene);
fun = @(x)betamle_23(x,ref(cursites),alt(cursites));
params(3,:)  = fmincon(fun,starts,[],[],[],[], ...
        [smallbnd smallbnd],[Inf Inf],[],options);

disp('init params:');
disp(params);

%========= run iterative procedure to estimate beta params

%============
gerrnimp = 0.001; % genotyping error probability for typed SNPs
gerrimp = 0.02; % genotyp error probability for imputed SNPs
serr = 0.001; % sequencing error probability
%============

% coompute per site the likelihood of the counts resulting from a
% genotyping error
gerrps = log(0.5)+ ...
         logsum_faster([ref*log(serr)+alt*log(1-serr), ...
                 alt*log(serr)+ref*log(1-serr)]);
gerrps(imputed(snp)) = log(gerrimp)+gerrps(imputed(snp));
gerrps(~imputed(snp)) = log(gerrnimp)+gerrps(~imputed(snp));

% obtain genotyping error probability per site (depends on typed/imputed)
nogerrps = zeros(size(ref));
nogerrps(imputed(snp)) = log(1-gerrimp);
nogerrps(~imputed(snp)) = log(1-gerrnimp);

siten = length(ref);
usesnp = true(siten,1); % snps to be used (=free of genotyping error)
sitenogerrps = zeros(siten,3); % (site x class) likelihoods assuming no genotyping error
siteps = zeros(siten,3); % total (site x class) probability
geneps = zeros(genn,3); % total (gene x class) probability
maxgene = max(gene);

ps = ones(1,3)/3; % inital class probabilities are even
itnum = 20;
for it = (1:itnum),
    
    tic;
    disp(['it: ' num2str(it)]);

    disp('computing zs ...');
    
    % compute site likelihood under each of the three models
    % assuming no genotyping error
    k = 1;
    sitenogerrps(:,k) =  betaln(ref+params(k,1),alt+params(k,2)) ...
                         - betaln(params(k,1),params(k,2));
    for k = [2 3],
        sitenogerrps(:,k) = log(0.5) ...
                            + logsum_faster([betaln(ref+params(k,1),alt+params(k,2)) ...
                                      - betaln(params(k,1),params(k,2)), ...
                                      betaln(alt+params(k,1),ref+params(k,2)) ...
                                      - betaln(params(k,1),params(k,2))]);  
    end
    
    % integrate genotyping error probabilities
    sitenogerrps = sitenogerrps + repmat(nogerrps,[1 3]);
    for k = (1:3),       
        siteps(:,k) = logsum_faster([gerrps sitenogerrps(:,k)]);
    end
    for k = (1:3),
         geneps(1:maxgene,k) = accumarray(gene,siteps(:,k));  
    end

    % integrate per-SNP likelihoods to per-gene probabilities
    for k = (1:3),
        geneps(:,k) = geneps(:,k) + log(ps(k));  
    end
    geneps = normrows_faster(geneps);
    geneps = geneps./repmat(sum(geneps,2),[1 3]);
    
    % draw class per gene
    tmp = mnrnd(1,geneps);
    zs = (tmp==1);
    [trash tmp] = sort(tmp,2,'descend');
    
    tmp = tmp(:,1); % tmp is a vector of class numbers per gene
    tmp = tmp(gene); % tmp is a vector of class numbers per site

    % compute per-site genotyping error probabilities
    tmp = sitenogerrps(sub2ind([siten 3],(1:siten)',tmp));
    tmp = normrows_faster([gerrps tmp]);
    tmp = tmp./repmat(sum(tmp,2),[1 2]);
    % draw genotyping error state per site
    usesnp = mnrnd(1,tmp);
    usesnp = (usesnp(:,2)==1);

    % estimate parameters given assigned classes and error states

    disp('optimizing ...');
    k=1;
    cursites = ismember(gene,find(zs(:,k))) & usesnp;
    fun = @(x)betamle_1(x,ref(cursites),alt(cursites));
    params(k,:) = fmincon(fun,params(k,:),[],[],[],[], ...
        [smallbnd smallbnd],[Inf Inf],[],options);    
    for k = [2 3],
        cursites = ismember(gene,find(zs(:,k))) & usesnp;
        fun = @(x)betamle_23(x,ref(cursites),alt(cursites));
        params(k,:) = fmincon(fun,params(k,:),[],[],[],[], ...
            [smallbnd smallbnd],[Inf Inf],[],options);
    end

    ps = sum(zs)/sum(sum(zs));

    disp('current parameters:');
    disp(params);
    
    toc;
end % endof foreach it

disp('final parameter values are:');
disp(params);

