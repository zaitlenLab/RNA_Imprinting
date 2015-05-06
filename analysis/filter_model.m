% This code runs the filtering and imprinting analysis steops on the
% Geuvadis LCLs dataset

% !! replace with the directory into which imprinting_code was extracted
refdir = '/Users/kthefrog/path_to_dir/imprinting_code';

%=====================
% data files and flags
%=====================

% directories of code files
addpath([refdir '/analysis']);
addpath([refdir '/common']);

% directory of data files
dprefix = [refdir '/data'];
auxprefix = [refdir '/aux_files'];
rprefix = [refdir '/results'];

% dataset used for inference
% datset for sample run:
datafile = 'GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.LCL.SAMPLE';
% !! replace with this file for the complete Geuvadis datset
%datafile = 'GD462.ASE.COV8.ANNOTPLUS.SINFO.txt.LCL';

% (individual,gene) pairs for which NMD event is predicted
% these will be omitted from the analysis
nmdfile = [auxprefix '/' 'GD_ind_gene_NMD.txt'];
% a list of types (as opposed to imputed) SNPs
typedfile = [auxprefix '/' 'GEUVADIS.PH1PH2_465.GTd.snps.txt'];
% a list of individual from whom phasing in not given
gdphasefile = [auxprefix '/' 'imputed_samples_phase_error.txt'];

%=================
% annotation files
%=================

% list of eqtl effects;
eqtlfile = [auxprefix '/' 'GD_eQTL_ase_effects.txt'];
[eqtl_gen t t t eqtl_tis t t eqtl_var t] = ...
    textread(eqtlfile,'%s%s%s%s%s%s%s%s%s','headerlines',1);

% list of known imprinted genes
% taken geneimprint.com and Morcos et al 2011
knownfile = [auxprefix '/' 'known_morcos_geneimprint_su.txt'];
knimp_gen = textread(knownfile,'%s');

% list of random monoallelic expression genes
% taken from Gimelbrant et al 2007
rmefile = [auxprefix '/' 'rme_gimelbrant_ensg_su.txt'];
rme_gen = textread(rmefile,'%s');

% conversion from ensg to name
ensg2namefile = [auxprefix '/' 'ensg2name.txt'];
[ensg gname] = textread(ensg2namefile,'%s%s','delimiter',',');

%==========
% READ DATA
%==========

disp(['reading data for ' datafile]);

%========= load counts and phase

load([dprefix '/' datafile '.mat'],'ref','n','ph');
ph(ph==-1) = 0;
ph = (ph==1);
snpn = size(n,2);

%========= load genes

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

[trash ind] = sortrows([chr(i2o_rs) pos(i2o_rs)],[1 2]);
u_rs = u_rs(ind);
typed = textread(typedfile,'%s');
imputed = ~ismember(u_rs,typed);
clear pos ind trash chr; 

%=======
% FILTER
%=======

disp('filtering');

remove = false(snpn,1);

%========= HWE FILTERING

load([dprefix '/' datafile '.hwe.mat'],'hwep');
hweremove = (0<=hwep & hwep<1e-3)';
remove = (remove|hweremove);

%========= FLIP FILTERING - PER SNP

thr2 = 1e-3;
tmp1 = 16<=n;
tmp1(:,remove) = false;
tmp1 = tmp1(:);
monoserr = 0.001;
% parameters used in the "flip test" for calling monoallelic sites
% identical to the beta parameters used for computing the per-gene
% likelihoods per each of the three expression classes
monoparams = [ ...
   45.5083   44.5147
    6.3187    5.6951
    0.6355    0.1493
];

siten = sum(sum(tmp1));
sitel = -Inf*ones(siten,3,2); % siten x modeln x allelen

% compute per site the likelihood under each of the three models
% (balanced, imbalanced, imprinted)

for model = 1,
    sitel(:,model,1) = betaln(monoparams(model,1)+ref(tmp1), ...
                    monoparams(model,2)+n(tmp1)-ref(tmp1)) ...
                    -betaln(monoparams(model,1),monoparams(model,2));
end

for model = [2 3],
    sitel(:,model,1) = betaln(monoparams(model,1)+ref(tmp1), ...
                    monoparams(model,2)+n(tmp1)-ref(tmp1)) ...
                    -betaln(monoparams(model,1),monoparams(model,2));
    sitel(:,model,2) = betaln(monoparams(model,1)+n(tmp1)-ref(tmp1), ...
                    monoparams(model,2)+ref(tmp1)) ...
                    -betaln(monoparams(model,1),monoparams(model,2));
end

sitel = max(sitel,[],3);
ldiff = sitel(:,3)-max(sitel(:,1:2),[],2);
% monoallelic sites are those for which max(like1lik2) < like3
mono = ldiff>0;
if (sum(mono)>0),
    f = ref(tmp1)./n(tmp1);
    f = f(:);
    tmpsnp = repmat(1:size(ref,2),[size(ref,1) 1]);
    whichsnp = tmpsnp(tmp1);
    whichsnp = whichsnp(:);
    clear tmpsnp;
    tmp3 = zeros(size(ref,2),2);
    % "flip test": use a binomial test to detect genes in which monoallelic
    % SNPs have a preference for ref or alt
    tmp3(1:max(whichsnp(mono)),1) = accumarray(whichsnp(mono),f(mono)<0.5);
    tmp3(1:max(whichsnp(mono)),2) = accumarray(whichsnp(mono),f(mono)>0.5);
    [trash tmp4] = binofit(tmp3(:,1),sum(tmp3,2),thr2);

    flipremove = tmp4(:,2)<0.5 | 0.5<tmp4(:,1);
    bothmono = all(tmp3>0,2);
    remove = (remove|flipremove);
else
    disp('no monos; skipping flip');
    bothmono = false(size(ref,2),1)
end
clear tmp1 whichsnp tmp3 tmp4 f sitel monoparams monoserr;

%===============

if (strcmp(u_gen{1},'')),
	remove = (remove|i2u_gen==1);
end

% clean removed SNPs
ref = ref(:,~remove);
n = n(:,~remove);
ph = ph(:,~remove);
i2u_gen = i2u_gen(~remove);
imputed = imputed(~remove);
u_rs = u_rs(~remove);
bothmono = bothmono(~remove);

% clean NMD sites
nmd = nmd(:,~remove);
n(nmd) = 0;
ref(nmd) = 0;
ph(nmd) = false;
clear nmd;

% keep only one SNP per gene in the small groups of individuals for whom no
% phasing information is available
gdphase = textread(gdphasefile,'%s');
ind2fix = find(ismember(u_id,gdphase));
for thisind = ind2fix',
    nonan = (n(thisind,:)>0);
    rlvgenes = unique(i2u_gen(nonan));
    for thisgen = rlvgenes',
        useful = find(nonan' & i2u_gen==thisgen);
        if (1<length(useful)),
             useful(randi(length(useful),1)) = [];
             n(thisind,useful) = 0;
             ref(thisind,useful) = 0;                 
             ph(thisind,useful) = false;                 
        end
    end
end

clear rs hwep hweremove i2o_gen i2o_rs i2o_id id ind2fix;

%===================
% MODEL COMPUTATIONS
%===================

disp('performing likelihood computations');

params = [ ...
   45.5083   44.5147
    6.3187    5.6951
    0.6355    0.1493
];

%============
gerrnimp = 0.001; % genotyping error rate in nonimputed SNPs
gerrimp = 0.05; % genotyping error rate in imputed SNPs
serr = 0.001; % sequencing error rate
perr = 0.2; % phasing error rate
%============

likes = ones(indn,3,genn); % likelihood of each class, per (ind x gene)
nonan = (0<n);
hascov = any(nonan,1);

totalindn = zeros(genn,1);
totalsnpn = zeros(genn,1);
snpthr = 3;
goodsnpn = zeros(genn,1);
indthr = 2;
goodindn = zeros(genn,1);
biasf = zeros(genn,4);
reff = zeros(genn,2);

snpperind = zeros(genn,indn); % number of informative SNPs per (ind x gene) 

hasboth = zeros(genn,1);
for i = (1:genn), % for every gene

    hasboth(i) = any(bothmono(i2u_gen==i));
    
    thissnp = (i2u_gen==i)' & hascov;
    thissnp = find(thissnp);
    thisind = any(nonan(:,thissnp),2);
    thisind = find(thisind);
    if (isempty(thisind)),
        continue;
    end
    thisnonan = nonan(thisind,thissnp);
    thisindn = length(thisind);
    thissnpn = length(thissnp);
        
    snpperind(i,thisind) = sum(thisnonan,2);
    
    thisref = ref(thisind,thissnp);
    thisref = thisref(thisnonan);
    thisref = thisref(:);
    thisalt = n(thisind,thissnp)-ref(thisind,thissnp);
    thisalt = thisalt(thisnonan);
    thisalt = thisalt(:);
    
    thisph = ph(thisind,thissnp);
    thisph = thisph(thisnonan);
    thisph = thisph(:);
    
    thisind = repmat(thisind,[1 thissnpn]);
    thisind = thisind(thisnonan);
    thisind = thisind(:);
    [u_ind i2o_ind i2u_ind] = unique(thisind);

    thissnp = repmat(thissnp,[thisindn 1]);
    thissnp = thissnp(thisnonan);
    thissnp = thissnp(:);
    [u_snp i2o_snp i2u_snp] = unique(thissnp); % just for useful stat; remove
    
    thisimp = imputed(thissnp);
       
    %=============================
    % collect some useful information
    totalindn(i) = thisindn;
    goodindn(i) = sum(indthr<=sum(thisnonan,2));
    totalsnpn(i) = thissnpn;
    goodsnpn(i) = sum(snpthr<=sum(thisnonan,1));
    %
    tmp = thisref./(thisref+thisalt);
    reff(i,1) = mean(tmp);
    reff(i,2) = var(tmp);
    tmp = max([tmp 1-tmp],[],2);
    biasf(i,1) = mean(tmp);
    biasf(i,2) = var(tmp);
    tmp1 = accumarray(i2u_ind,thisref);
    tmp2 = accumarray(i2u_ind,thisalt);
    tmp = tmp1./(tmp1+tmp2);
    tmp = max([tmp 1-tmp],[],2);
    biasf(i,3) = mean(tmp);
    biasf(i,4) = var(tmp);
    %=============================
    
    % compute the likelihood of every site being a genotyping error
    thisgerr = zeros(size(thisref));
    thisgerr(thisimp) = gerrimp;
    thisgerr(~thisimp) = gerrnimp;
    thisgerrps = log(thisgerr)+log(0.5)+logsum_faster([thisref*log(serr)+thisalt*log(1-serr), ...
                                                       thisref*log(1-serr)+thisalt*log(serr)]);
    tmpz = zeros(thisindn,3);
    
    % compute model1 (balanced) likelihood
    model = 1;
    tmpbal = betaln(params(model,1)+thisref, ...
                    params(model,2)+thisalt) ...
            -betaln(params(model,1),params(model,2));
    tmpbal = logsum_faster([thisgerrps log(1-thisgerr)+tmpbal]);
    
    tmpz(:,model) = accumarray(i2u_ind,tmpbal); 
    
    thisrefbu = thisref;
    thisref(thisph) = thisalt(thisph);
    thisalt(thisph) = thisrefbu(thisph);
    
    % compute model2 (imbalanced) and model3 (imprinted)
    % likelihoods
    for model = [2 3],
            % compute per site the model likelihoods for the two cases: 
            % gene copy #1 is over-expressed
            % gene copy #2 is over-expressed
			c1exp = betaln(params(model,1)+thisref,params(model,2)+thisalt) ...
				   -betaln(params(model,1),params(model,2));
			c1sil = betaln(params(model,1)+thisalt,params(model,2)+thisref) ...
				   -betaln(params(model,1),params(model,2));
        % join the site likelihood into a single lieklihood per gene copy
        % while allowing for phasing errors
        c1expwph = logsum_faster([log(1-perr)+c1exp log(perr)+c1sil]);
        c1silwph = logsum_faster([log(1-perr)+c1sil log(perr)+c1exp]);
        c1expwphwg = logsum_faster([thisgerrps log(1-thisgerr)+c1expwph]);
        c1silwphwg = logsum_faster([thisgerrps log(1-thisgerr)+c1silwph]);
        tmp = zeros(thisindn,2);
        tmp(:,1) = accumarray(i2u_ind,c1expwphwg);
        tmp(:,2) = accumarray(i2u_ind,c1silwphwg);
        tmpz(:,model) = log(0.5)+logsum_faster(tmp);
    end

    likes(u_ind,:,i) = tmpz;
    
end

%========= COMPUTE MIXTURE LIKELIHOODS

disp('computing imprinting statistics');

%==========
minsnp = 1;
%==========

allimpl = zeros(genn,1);
noimplem = zeros(genn,1);
someimplem = zeros(genn,1);
optxem = zeros(genn,3);
for i = (1:genn),
    
    thisind = minsnp<=snpperind(i,:);
    if (sum(thisind)==0)
        continue;
    end
        
    tmp = likes(thisind,:,i);
    
    allimpl(i) = sum(tmp(:,3));

    [trash noimplem(i)] = em12(tmp(:,1:2),ones(1,2)/2);
    
    [optxem(i,:) someimplem(i)] = em12(tmp,ones(1,3)/3);
    
end

impglr = someimplem-noimplem;
hetlr = someimplem-max([allimpl noimplem],[],2);
implr = allimpl-noimplem;

impglr(impglr<0) = 0;
hetlr(hetlr<0) = 0;

%=========  COMPUTE PER-GENE AND PER-GENE, PER-INDIVIDUAL Z'S

count = zeros(genn,3,2); % is computed using snpn >= minsnp
indcount = zeros(genn,2);
zs = -1*ones(indn,3,genn); % is computed for snpn >= 1
gzs = -1*ones(genn,3); % is computed using snpn >= 1

%==========
ps = [0.8 0.19 0.01]; % prior class probabilities
impthr = 0.7; % probability threshold for an individual to be classified
%==========

for i = (1:genn),

    thisind = 1<=snpperind(i,:);
    if (sum(thisind)==0)
        continue;
    end
    zs(thisind,:,i) = ...
        normrows_faster(likes(thisind,:,i)+repmat(log(ps),[sum(thisind) 1]));
    gzs(i,:) = normrows_faster(sum([log(ps); likes(thisind,:,i)]));
    
    for minsnp = [1 2],

        indcount(i,minsnp) = sum(minsnp<=snpperind(i,:));
        count(i,:,minsnp) = sum(zs(minsnp<=snpperind(i,:),:,i)>=impthr);

    end

end

%=========================    
% PRINT GENE LIST AND FIGS
%=========================    

disp('outputting results');

%========= OBTAIN GENE INFO

eqtl_str = cell(genn,2);
eqtl_str(:) = {'-'};
for i = (1:genn),
    thisind = find(strcmp(eqtl_gen,u_gen(i)),1,'first');
    if (~isempty(thisind)),
       eqtl_str{i,1} = eqtl_tis{thisind};
       eqtl_str{i,2} = eqtl_var{thisind};
    end
end

knimp = ismember(u_gen,knimp_gen);
knimp = knimp+1;

rme = ismember(u_gen,rme_gen);
rme = rme+1;

hasboth = hasboth+1;

name_str = cell(1,genn);
name_str(:) = {'-'};
for i = (1:genn),
    thisind = find(strcmp(ensg,u_gen(i)),1,'first');
    if (~isempty(thisind)),
       name_str{i} = gname{thisind};
    end
end


%========= WRITE RESULTS, SORTED ACCORDING TO IMPGLR

[trash ind] = sort(impglr,'descend');
ind = ind(ismember(ind,find(sum(snpperind,2)>0)));
yesno = ['-','+'];

fid = fopen([rprefix '/' datafile '.impglr.res.txt'],'w');

fprintf(fid,['%-4s\t' ... % rank
             '%-20s' ...    % ensg
             '%-20s' ...    % name             
             '%-15s' ...    % impglr    
             '%-15s' ...    % implr                 
             '%-15s' ...    % hetlr 
             '%-10s' ...    % opz1p              
             '%-10s' ...    % opz2p              
             '%-10s' ...    % opz3p              
             '%-5s' ...    % bal inds >= 2
             '%-5s' ...    % imbal inds >= 2
             '%-5s' ...    % imprinted inds >= 2
             '%-8s' ...    % total inds >= 2
             '%-5s' ...    % bal inds >= 1
             '%-5s' ...    % imbal inds >= 1
             '%-5s' ...    % imprinted inds >= 1
             '%-8s' ...    % total inds >= 1
             '%-6s' ...    % gz1
             '%-6s' ...    % gz2
             '%-6s' ...    % gz3
             '%-7s' ...    % hasboth
             '%-7s%-7s\t\t\t%-80s%-80s' ...  % known, rme, eqtal tis, eqtl val
             '%-10s' ...    % indn
             '%-10s' ...    % snpn
             '%-10s' ...    % indn2
             '%-10s' ...    % snpn3
             '%-10s' ...    % biasf(1)
             '%-10s' ...    % biasf(2)
             '%-10s' ...    % biasf(3)
             '%-10s' ...    % biasf(4)
             '%-10s' ...    % reff(1)
             '%-10s' ...    % reff(2)
             '\n'], ... % eqtal pval, eqtl tis, known
            'rank','ensg','name', ...
            'impglr','implr','hetlr', ...
            'optz1p','optz2p','optz3p', ...
            'bal2','imb2','imp2','total2', ...
            'bal1','imb1','imp1','total1', ...
            'gz1','gz2','gz3', ...
            'both', ...
            'known','RME','eQTLt','eQTLv', ...
            'indn','snpn','indn2','snpn3', ...
            'max_m','max_v','maxph_m','maxph_v', ...
            'ref_m','ref_v');
        
for i = (1:length(ind)),

fprintf(fid,['%4d\t' ... % rank
             '%-20s' ...    % ensg
             '%-20s' ...    % name
             '%-15.3f' ...    % impglr
             '%-15.3f' ...    % impglr
             '%-15.3f' ...    % hetlr             
             '%-10.2f' ...    % optz1p
             '%-10.2f' ...    % optz2p
             '%-10.2f' ...    % optz3p
             '%-5d' ...    % bal2 count
             '%-5d' ...    % imb2 count
             '%-5d' ...    % imp2 count
             '%-8d' ...    % ind2 count
             '%-5d' ...    % bal1 count
             '%-5d' ...    % imb1 count
             '%-5d' ...    % imp1 count
             '%-8d' ...    % ind1 count
             '%-6.2f' ...    % gz1
             '%-6.2f' ...    % gz2
             '%-6.2f' ...    % gz3
             '%-7s' ...    % hasboth
             '%-7s%-7s\t\t\t%-80s%-80s' ... % known, rme, eqtal tis, eqtl val
             '%-10d' ...    % indn
             '%-10d' ...    % snpn
             '%-10d' ...    % indn2
             '%-10d' ...    % snpn3
             '%-10.3f' ...    % biasf(1)
             '%-10.3f' ...    % biasf(2)
             '%-10.3f' ...    % biasf(3)
             '%-10.3f' ...    % biasf(4)
             '%-10.3f' ...    % reff(1)
             '%-10.3f' ...    % reff(2)
             '\n'], ... 
            i, ... 
            u_gen{ind(i)}, ...
            name_str{ind(i)}, ...
            impglr(ind(i)), ...            
            implr(ind(i)), ...            
            hetlr(ind(i)), ...
            optxem(ind(i),1), ...
            optxem(ind(i),2), ...
            optxem(ind(i),3), ...
            count(ind(i),1,2), ...
            count(ind(i),2,2), ...
            count(ind(i),3,2), ...
            indcount(ind(i),2), ...
            count(ind(i),1,1), ...
            count(ind(i),2,1), ...
            count(ind(i),3,1), ...
            indcount(ind(i),1), ...
            gzs(ind(i),1), ...
            gzs(ind(i),2), ...
            gzs(ind(i),3), ...
            yesno(hasboth(ind(i))), ...
            yesno(knimp(ind(i))),yesno(rme(ind(i))),eqtl_str{ind(i),1},eqtl_str{ind(i),2}, ...
            totalindn(ind(i)), ...
            totalsnpn(ind(i)), ...
            goodindn(ind(i)), ...
            goodsnpn(ind(i)), ...
            biasf(ind(i),1), ...    
            biasf(ind(i),2), ...    
            biasf(ind(i),3), ...    
            biasf(ind(i),4), ...    
            reff(ind(i),1), ...    
            reff(ind(i),2));

end

fclose(fid);

% print scatter plots for (at most) 20 top genes
counter = 0;
i = 1;
classnames = {'BAL','IMB','IMP'};
while (counter<20 && i<=length(ind)),

    if (hasboth(ind(i))==1),
        i = i+1;
        continue;
    end
    
    thissnp = (i2u_gen==ind(i))' & hascov;
    thissnp = find(thissnp);
    thisind = any(nonan(:,thissnp),2);
    thisind = find(thisind);
    if (isempty(thisind)), % unnec
        i = i+1;
        continue;
    end
    thisnonan = nonan(thisind,thissnp);
    thisindn = length(thisind);
    thissnpn = length(thissnp);

    thisref = ref(thisind,thissnp);
    thisref = thisref(thisnonan);
    thisref = thisref(:);
    thisalt = n(thisind,thissnp)-ref(thisind,thissnp);
    thisalt = thisalt(thisnonan);
    thisalt = thisalt(:);

    thisind = repmat(thisind,[1 thissnpn]);
    thisind = thisind(thisnonan);
    thisind = thisind(:);
    [u_ind i2o_ind i2u_ind] = unique(thisind);

    for j = (1:3),
        clf;
        hold on;
        title({name_str{ind(i)},[' KNOWN' yesno(knimp(ind(i))) ' RME' yesno(rme(ind(i))) '  ']});
        scatter(thisref,thisalt,60,zs(u_ind(i2u_ind),j,ind(i)),'filled');
        xlim([-1 max(thisref)+1]);
        ylim([-1 max(thisalt)+1]);
        xlabel('ref counts');
        ylabel('alt counts');      
        colorbar;
        axis square;
        saveas(gcf,[rprefix '/' datafile '.impglr.res.hit.' num2str(i) '.' name_str{ind(i)} '.' classnames{j} '.jpg']);
    end

    i = i+1;
    counter = counter+1;

end    

