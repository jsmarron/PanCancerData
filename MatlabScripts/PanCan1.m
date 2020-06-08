disp('Running MATLAB script file PanCan1.m') ;
%
%    For analysis of TCGA PanCan data set from Chris Fan
%    12 Cancer Types
%
%    The main input file is:    kviao.pancan12tu.txt
%        Note 1st three lines have text:
%            - Long names
%            - Cancer Type
%            - Platform
%        Remaining lines 4 - 12485 are Gene Expression
%


ipart = 16 ;    %  0 - Read data, check inputs, store as .mat file
               %  1 - Check Missings
               %  2 - Check Data normalizations (only genes w/ no missings)
               %  3 - Make PCA Scatterplots (only genes w/ no missings)
               %  4 - Explore Mean MargDist plots
               %  5 - Explore Subset COAD, READ, UCEC
               %  6 - Classify each type versus the rest
               %  10 - Do better handling of missings, and store as .mat file
               %  11 - PCA Scatterplots
               %  12 - Explore Mean MargDist plots
               %  13 - Explore Subset COAD, READ, UCEC
               %  14 - Classify each type versus the rest
               %  15 - MargOverlap, AUC-ROC, each versus the rest 
               %  16 - Explore LUAD vs. LUSC
               %  20 - Pairwise classifications



savefilestr = 'PanCan1.mat' ;
savefileAstr = 'PanCan1A.mat' ;



if ipart == 0 ;    %  Check inputs

  %  Read in Main Data Values
  %
  disp(' ') ;
  disp('    Start Main Data Load') ;
  tic ;
  infilestr = 'kviao.pancan12tu.txt' ;
  instruct = importdata(infilestr,'\t',3) ;
      %  treat first 19 lines as headers
  disp('    Finished Main Data Load') ;
  disp(['        took ' num2str(toc / 60) ' minutes']) ;

  mdata = instruct.data ;
  catext = instruct.textdata ;

  %  Get gene names
  %
  caGeneName = catext(4:end,1) ;
      %  Names of variables (with numerical data) from first column of text input
  disp(' ') ;
  disp(' ') ;
  disp('    Check Text Inputs') ;
  disp(' ') ;
  disp('    Check First Gene Name is: ABP1|26') ;
  caGeneName{1}
  disp(' ') ;
  disp('    Check Last Gene Name is: TARDBP|23435') ;
  caGeneName{end}

  %  Get CLIDs
  %
  caname = textscan(catext{1,1},'%s') ;
  ninput = size(caname{1},1) ;
  caCLIDlong = caname{1}([4:5:ninput]) ;
      %  Long version of CLIDs
      %  Need to trim to 12 characters
  disp('        Number of cases = ') ;
  n = size(caCLIDlong,1) 
  for i = 1:n ;
    caCLID{i} = caCLIDlong{i}(1:12) ;
  end ;
  disp(' ') ;
  disp('    Check CLID Inputs') ;
  disp(' ') ;
  disp('    Check First CLID is: TCGA-DK-A3WX') ;
  caCLID{1}
  disp(' ') ;
  disp('    Check Last CLID is: TCGA-D8-A1J9') ;
  caCLID{end}

  %  Get Cancer Types
  %
  caCanType = caname{1}([2:5:ninput]) ;
  disp(' ') ;
  disp('    Check number of Cancer Types is right (this is 0):') ;
  n - size(caCanType,1) 
  disp(' ') ;
  disp('    Check CanType Inputs') ;
  disp(' ') ;
  disp('    Check First CanType is: BLCA') ;
  caCanType{1}
  disp(' ') ;
  disp('    Check Last CanType is: BRCA') ;
  caCanType{end}

  %  Get Platforms
  %
%  caplat = textscan(catext{3,1},'%s','Delimiter','\t') ;
%  caplat = caplat{1} ;
  disp(' ') ;
  disp('    Check Platform Inputs') ;
  disp(' ') ;
  disp('    Check Header is: Platform') ;
  catext{3,1} ;
  caPlatform = {} ;
  for i = 1:n ;
    caPlatform{i} = catext{3,i+1} ;
  end ;
  disp(' ') ;
  disp('    Check number of forms is right (this is 0):') ;
  n - size(caPlatform,1) 
  disp(' ') ;
  disp('    Check First Platform is: HiSeq') ;
  caPlatform{1}
  disp(' ') ;
  disp('    Check Last Platform is: HiSeq') ;
  caPlatform{end}


  disp(' ') ;
  disp(' ') ;
  disp('    Check Numerical Data Inputs') ;
  disp('    Check number of inputs is right (this is 0):') ;
  n - size(mdata,2) 
  disp(' ') ;
  disp('    Check Numerical Data') ;
  disp(' ') ;
  disp('    Check mdata(1,1) is 9.6:') ;
  mdata(1,1)
  disp('    Check mdata(1,end) is 4.24:') ;
  mdata(1,end)
  disp('    Check mdata(end,1) is 11.92:') ;
  mdata(end,1)
  disp('    Check mdata(end,end) is 11.74:') ;
  mdata(end,end)
  disp(' ') ;


  disp(' ') ;
  disp('    Start .mat file save') ;
  save(savefilestr,'mdata','caGeneName','caCLID','caCanType','caPlatform') ;
  disp('    Finished .mat file save') ;
  disp(' ') ;


else ;    %  Work with .mat version of data

  disp(' ') ;
  disp('    Start .mat file read') ;
  load(savefilestr) ;
  %    Loads Variables:
  %        mdata:         Matrix of Gene Expression
  %        caGeneName:    Cell Array of Gene Names
  %        caCLID         Cell Array of Case identifiers
  %        caCanType      Cell Array of Cancer Types
  %        caPlatform     Cell Array of Platforms

  n = size(mdata,2) ;
  d = size(mdata,1) ;

  mcolorbase = RainbowColorsQY(12) ;
  caGeneNamebase = {'BLCA' 'COAD' 'GBM' 'KIRC' 'LAML' 'LUAD' ...
                        'LUSC' 'OV' 'READ' 'HNSC' 'UCEC' 'BRCA'} ;

  %  Sort out Subtypes
  %
  vBLCAflag = zeros(1,n) ;
  vCOADflag = zeros(1,n) ;
  vGBMflag = zeros(1,n) ;
  vKIRCflag = zeros(1,n) ;
  vLAMLflag = zeros(1,n) ;
  vLUADflag = zeros(1,n) ;
  vLUSCflag = zeros(1,n) ;
  vOVflag = zeros(1,n) ;
  vREADflag = zeros(1,n) ;
  vHNSCflag = zeros(1,n) ;
  vUCECflag = zeros(1,n) ;
  vBRCAflag = zeros(1,n) ;
  mcolor = zeros(n,3) ;
  for i = 1:n ;
    if strcmp(caCanType{i},'BLCA') ;
      vBLCAflag(i) = 1 ;
      mcolor(i,:) = mcolorbase(1,:) ;
    elseif strcmp(caCanType{i},'COAD') ;
      vCOADflag(i) = 1 ;
      mcolor(i,:) = mcolorbase(2,:) ;
    elseif strcmp(caCanType{i},'GBM') ;
      vGBMflag(i) = 1 ;
      mcolor(i,:) = mcolorbase(3,:) ;
    elseif strcmp(caCanType{i},'KIRC') ;
      vKIRCflag(i) = 1 ;
      mcolor(i,:) = mcolorbase(4,:) ;
    elseif strcmp(caCanType{i},'LAML') ;
      vLAMLflag(i) = 1 ;
      mcolor(i,:) = mcolorbase(5,:) ;
    elseif strcmp(caCanType{i},'LUAD') ;
      vLUADflag(i) = 1 ;
      mcolor(i,:) = mcolorbase(6,:) ;
    elseif strcmp(caCanType{i},'LUSC') ;
      vLUSCflag(i) = 1 ;
      mcolor(i,:) = mcolorbase(7,:) ;
    elseif strcmp(caCanType{i},'OV') ;
      vOVflag(i) = 1 ;
      mcolor(i,:) = mcolorbase(8,:) ;
    elseif strcmp(caCanType{i},'READ') ;
      vREADflag(i) = 1 ;
      mcolor(i,:) = mcolorbase(9,:) ;
    elseif strcmp(caCanType{i},'HNSC') ;
      vHNSCflag(i) = 1 ;
      mcolor(i,:) = mcolorbase(10,:) ;
    elseif strcmp(caCanType{i},'UCEC') ;
      vUCECflag(i) = 1 ;
      mcolor(i,:) = mcolorbase(11,:) ;
    elseif strcmp(caCanType{i},'BRCA') ;
      vBRCAflag(i) = 1 ;
      mcolor(i,:) = mcolorbase(12,:) ;
    else ;
      disp('    Error from PanCan1.m') ;
      disp(['    Unknown Cancer Type:  ' caCanType{i}]) ;
      pauseSM ;
    end ;
  end ;
  vBLCAflag = logical(vBLCAflag) ;
  vCOADflag = logical(vCOADflag) ;
  vGBMflag = logical(vGBMflag) ;
  vKIRCflag = logical(vKIRCflag) ;
  vLAMLflag = logical(vLAMLflag) ;
  vLUADflag = logical(vLUADflag) ;
  vLUSCflag = logical(vLUSCflag) ;
  vOVflag = logical(vOVflag) ;
  vREADflag = logical(vREADflag) ;
  vHNSCflag = logical(vHNSCflag) ;
  vUCECflag = logical(vUCECflag) ;
  vBRCAflag = logical(vBRCAflag) ;
  mCTflag = [vBLCAflag; ...
           vCOADflag; ...
           vGBMflag; ...
           vKIRCflag; ...
           vLAMLflag; ...
           vLUADflag; ...
           vLUSCflag; ...
           vOVflag; ...
           vREADflag; ...
           vHNSCflag; ...
           vUCECflag; ...
           vBRCAflag] ;

  disp(' ') ;
  disp(['    Number of BLCA cases is:  ' num2str(sum(vBLCAflag))]) ;
  disp(['    Number of COAD cases is:  ' num2str(sum(vCOADflag))]) ;
  disp(['    Number of GBM cases is:  ' num2str(sum(vGBMflag))]) ;
  disp(['    Number of KIRC cases is:  ' num2str(sum(vKIRCflag))]) ;
  disp(['    Number of LAML cases is:  ' num2str(sum(vLAMLflag))]) ;
  disp(['    Number of LUAD cases is:  ' num2str(sum(vLUADflag))]) ;
  disp(['    Number of LUSC cases is:  ' num2str(sum(vLUSCflag))]) ;
  disp(['    Number of OV cases is:  ' num2str(sum(vOVflag))]) ;
  disp(['    Number of READ cases is:  ' num2str(sum(vREADflag))]) ;
  disp(['    Number of HNSC cases is:  ' num2str(sum(vHNSCflag))]) ;
  disp(['    Number of UCEC cases is:  ' num2str(sum(vUCECflag))]) ;
  disp(['    Number of BRCA cases is:  ' num2str(sum(vBRCAflag))]) ;




  if ipart == 1 ;    %  Check Missings

    mflagmiss = isnan(mdata) ;
    vnmissgene = sum(mflagmiss,2) ;
    vnmisscase = sum(mflagmiss,1) ;
    disp(' ') ;
    disp(['    # genes w/ missings = ' num2str(sum(vnmissgene > 0)) ...
                   '  out of ' num2str(d)]) ;
    disp(['    # cases w/ missings = ' num2str(sum(vnmisscase > 0)) ...
                   '  out of ' num2str(n)]) ;
    disp(['       Total # missings = ' num2str(sum(sum(mflagmiss,2))) ...
                   '  out of ' num2str(n * d)]) ;

    %  Plot # missings 
    %
    figure(1) ;
    clf ;
    plot((1:d)',vnmissgene,'k') ;
    title(['Explore missings per Gene, out of ' ...
                    num2str(n) ' cases']) ;
    xlabel('Gene Index') ;
    ylabel('# missing cases') ;
    orient landscape ;
    savestr = ['PanCan1ip' num2str(ipart) 'PlotNumMissing4Genes'] ;
    print('-dpsc2',savestr) ;

    figure(2) ;
    clf ;
    plot((1:n)',vnmisscase,'k') ;
    title(['Explore missings per Case, out of ' ...
                    num2str(d) ' genes']) ;
    xlabel('Case Index') ;
    ylabel('# missing genes') ;
    orient landscape ;
    savestr = ['PanCan1ip' num2str(ipart) 'PlotNumMissing4Cases'] ;
    print('-dpsc2',savestr) ;

    %  Above plots indicate 2 cases with many missing
    %
    %  First identify those
    %
    [vnmisscasesort,vicsort] = sort(vnmisscase,'descend') ;
    disp(' ') ;
    disp('    Biggest case counts are:') ;
    vnmisscasesort(1:5) 
    disp('    Biggest case indices are:') ;
    vicsort(1:5) 

    %  Create data matrix without 2 biggest cases
    %
    vflagnm1 = ones(1,n) ;
    vflagnm1(vicsort(1)) = 0 ;
    vflagnm1(vicsort(2)) = 0 ;
    vflagnm1 = logical(vflagnm1) ;
    mdatanm1 = mdata(:,vflagnm1) ;

    nnm1 = size(mdatanm1,2) ;
    mflagmissnm1 = isnan(mdatanm1) ;
    vnmissgenenm1 = sum(mflagmissnm1,2) ;
    disp(['    # genes w/ missings = ' num2str(sum(vnmissgenenm1 > 0)) ...
                   '  out of ' num2str(d)]) ;

    figure(3) ;
    clf ;
    plot((1:d)',vnmissgenenm1,'k') ;
    title(['Explore missings per Gene, out of ' ...
                    num2str(n) ' cases, without 2 cases']) ;
    xlabel('Gene Index') ;
    ylabel('# missing cases') ;
    orient landscape ;
    savestr = ['PanCan1ip' num2str(ipart) 'PlotNumMissing4Genes-nm1'] ;
    print('-dpsc2',savestr) ;

    %  Study distribution of missings per gene
    %
    [vnmissgenesort,vigsort] = sort(vnmissgenenm1,'descend') ;
    vindgsort = (1:d)' ;
    vindgsort = vindgsort(vigsort) ;

    figure(4) ;
    clf ;
    vflagn0 = vnmissgenesort > 0 ;
    dn0 = sum(vflagn0) ;
    plot((1:dn0)',vnmissgenesort(vflagn0),'k') ;
    title(['Explore missings per Gene, only missings,' ...
                   ' without 2 cases']) ;
    xlabel('Sorted Gene Ordering') ;
    ylabel('# missing cases') ;
    orient landscape ;
    savestr = ['PanCan1ip' num2str(ipart) 'PlotSortedNumMissing4Genes-nm1'] ;
    print('-dpsc2',savestr) ;

    %  Study Gene Non-Missing Standard Deviations
    %
    vgenesd = zeros(d,1) ;
    for i = 1:d ;
      vdata = mdata(i,:) ;
      vfmiss = isnan(vdata) ;
      sd = std(vdata(~vfmiss)) ;
      vgenesd(i) = sd ;
    end ;
    thresh = vnmissgenesort(2000) ;
    vflagbig = (vnmissgenenm1 > thresh) ;

    figure(5) ;
    clf ;
    plot(vnmissgenenm1(vindgsort(1:2000)),vgenesd(vindgsort(1:2000)),'bo') ;
    title(['Relationship between Gene SD and # missings']) ;
    xlabel('# missing cases') ;
    ylabel('Gene Standard Deviation (non-missings)') ;
    orient landscape ;
    savestr = ['PanCan1ip' num2str(ipart) 'MissingsVs SD'] ;
    print('-dpsc2',savestr) ;

    figure(6) ;
    clf ;
    dolcolor = ones(d,1) * [0.5 0.5 0.5] ;
    for i = 1:d ;
      if vflagbig(i) ;
        dolcolor(i,:) = [1 0 0] ;
      end ;
    end ;
    paramstruct = struct('ndatovlay',2, ...
                         'dolcolor',dolcolor, ...
                         'ylabelstr','Standard Deviation of Gene Expression') ;
    kdeSM(vgenesd,paramstruct) ;
    orient landscape ;
    savestr = ['PanCan1ip' num2str(ipart) 'KdeAllGeneSDs'] ;
    print('-dpsc2',savestr) ;



  elseif ipart == 2 ;    %  Check Data normalizations

    %  Work only with genes with no missings
    %
    mflagmiss = isnan(mdata) ;
    vnmissgene = sum(mflagmiss,2) ;
    vflagmissgene = (vnmissgene > 0) ;
        %  one for genes with a missing
    dnm = sum(~vflagmissgene) ;

    mdatanm = mdata(~vflagmissgene,:) ;
    caGeneNamenm = caGeneName(~vflagmissgene) ;

    %  MargDistPlot - Mean
    %
    figure(1) ;
    clf ;
    varnamecellstr = {caGeneNamenm} ;
    savestr = ['PanCan1ip' num2str(ipart) 'AllMargDistMean'] ;
    paramstruct = struct('istat',1, ...
                         'varnamecellstr',varnamecellstr, ...
                         'icolor',mcolor, ...
                         'isubpopkde',1, ...
                         'datovlaymax',0.9, ...
                         'datovlaymin',0.4, ...
                         'savestr',savestr) ;
    MargDistPlotSM(mdatanm,paramstruct) ;

    %  MargDistPlot - SD
    %
    figure(2) ;
    clf ;
    savestr = ['PanCan1ip' num2str(ipart) 'AllMargDistSD'] ;
    paramstruct = struct('istat',2, ...
                         'varnamecellstr',varnamecellstr, ...
                         'icolor',mcolor, ...
                         'isubpopkde',1, ...
                         'datovlaymax',0.9, ...
                         'datovlaymin',0.4, ...
                         'savestr',savestr) ;
    MargDistPlotSM(mdatanm,paramstruct) ;

    %  MargDistPlot - Skewness
    %
    figure(3) ;
    clf ;
    savestr = ['PanCan1ip' num2str(ipart) 'AllMargDistSkewness'] ;
    paramstruct = struct('istat',3, ...
                         'varnamecellstr',varnamecellstr, ...
                         'viplot',[(1:7), ((dnm - 7):dnm)], ...
                         'icolor',mcolor, ...
                         'isubpopkde',1, ...
                         'datovlaymax',0.9, ...
                         'datovlaymin',0.4, ...
                         'savestr',savestr) ;
    MargDistPlotSM(mdatanm,paramstruct) ;

    %  MargDistPlot - Kurtosis
    %
    figure(4) ;
    clf ;
    savestr = ['PanCan1ip' num2str(ipart) 'AllMargDistKurtosis'] ;
    paramstruct = struct('istat',4, ...
                         'varnamecellstr',varnamecellstr, ...
                         'viplot',((dnm - 14):dnm), ...
                         'icolor',mcolor, ...
                         'isubpopkde',1, ...
                         'datovlaymax',0.9, ...
                         'datovlaymin',0.4, ...
                         'savestr',savestr) ;
    MargDistPlotSM(mdatanm,paramstruct) ;


  elseif ipart == 3 ;    %  Make PCA Scatterplots

    %  Make Color Plot with Cancer Type Colors
    %
    figure(1) ;
    clf ;
    axis([0 1 0 1]) ;
    nct = 12 ;
    for ict = 1:nct ;
      text(0.3,(nct - ict + 1) / (nct + 1), ...
           caGeneNamebase{ict}, ...
           'Color',mcolorbase(ict,:), ...
           'FontSize',24) ;
      text(0.6,(nct - ict + 1) / (nct + 1), ...
           num2str(sum(mCTflag(ict,:))), ...
           'Color',mcolorbase(ict,:), ...
           'FontSize',24) ;
    end ;
    orient landscape ;
    savestr = ['PanCan1ip' num2str(ipart) 'CancerTypeColors'] ;
    print('-dpsc2',savestr) ;


    %  Work only with genes with no missings
    %
    mflagmiss = isnan(mdata) ;
    vnmissgene = sum(mflagmiss,2) ;
    vflagmissgene = (vnmissgene > 0) ;
        %  one for genes with a missing
    dnm = sum(~vflagmissgene) ;

    mdatanm = mdata(~vflagmissgene,:) ;

    %  PCA scatterplot, all together
    %
    figure(2) ;
    clf ;
    titlecellstr = {{'PanCan Gene Expression, 2013' ['n = ' num2str(n)] ...
                      ['# No-missings genes = ' num2str(dnm)] ...
                      ['out of total' num2str(d)]}} ;
    labelcellstr = {{'PC 1 Scores' 'PC 2 Scores' 'PC 3 Scores' 'PC 4 Scores'}} ;
    savestr = ['PanCan1ip' num2str(ipart) 'AllDataPCA-NonMissGenes'] ;
    paramstruct = struct('npcadiradd',4, ...
                         'icolor',mcolor, ...
                         'isubpopkde',1, ...
                         'ibelowdiag',0, ...
                         'titlecellstr',titlecellstr, ...
                         'labelcellstr',labelcellstr, ...
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    scatplotSM(mdatanm,[],paramstruct) ;

    %  Make some zoomed in plots
    %
    paramstruct = struct('npc',4,...
                         'iscreenwrite',1,...
                         'viout',[0 0 0 0 1]) ;
    outstruct = pcaSM(mdatanm,paramstruct) ;
    mpc = getfield(outstruct,'mpc') ;
        %  4 x n matrix of PC scores

    figure(3) ;
    clf ;
    savestr = ['PanCan1ip' num2str(ipart) 'AllData-NonMissGenes-PC1v2'] ;
    paramstruct = struct('icolor',mcolor, ...
                         'titlestr','Zoomed in PanCan All', ...
                         'xlabelstr','PC 2 Scores', ...
                         'ylabelstr','PC 1 Scores', ...
                         'ifigure',3, ...
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    projplot2SM(mpc,[[0; 1; 0; 0],[1; 0; 0; 0]],paramstruct) ;

    figure(4) ;
    clf ;
    savestr = ['PanCan1ip' num2str(ipart) 'AllData-NonMissGenes-PC2v3'] ;
    paramstruct = struct('icolor',mcolor, ...
                         'titlestr','Zoomed in PanCan All', ...
                         'xlabelstr','PC 3 Scores', ...
                         'ylabelstr','PC 2 Scores', ...
                         'ifigure',4, ...
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    projplot2SM(mpc,[[0; 0; 1; 0],[0; 1; 0; 0]],paramstruct) ;

    figure(5) ;
    clf ;
    savestr = ['PanCan1ip' num2str(ipart) 'AllData-NonMissGenes-PC3v4'] ;
    paramstruct = struct('icolor',mcolor, ...
                         'titlestr','Zoomed in PanCan All', ...
                         'xlabelstr','PC 4 Scores', ...
                         'ylabelstr','PC 3 Scores', ...
                         'ifigure',5, ...
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    projplot2SM(mpc,[[0; 0; 0; 1],[0; 0; 1; 0]],paramstruct) ;


  elseif ipart == 4 ;    %   Individual Mean MargDist plots & PCAs

    %  First do median Recentering of each gene
    %
    mdatamrc = [] ;
    for i = 1:d ;
      vdata = mdata(i,:) ;
      vflagmiss = isnan(vdata) ;
      med = median(vdata(~vflagmiss)) ;
      vdatamrc = vdata - med ;
      mdatamrc = [mdatamrc; vdatamrc] ;
    end ;

    disp('Finished Median Recentering') ;

    for ict = 1:12 ;

      strct = caGeneNamebase{ict} ;
      mdatamrcct = mdatamrc(:,mCTflag(ict,:)) ;
      mcolorct = mcolor(mCTflag(ict,:),:) ;
      nct = sum(mCTflag(ict,:)) ;

      disp(['    Working on Cancer Type ' strct]) ;

      mflagmiss = isnan(mdatamrcct) ;
      vnmissgene = sum(mflagmiss,2) ;
      vflagmissgene = (vnmissgene > 0) ;
          %  one for genes with a missing
      dnm = sum(~vflagmissgene) ;

      mdatamrcctnm = mdatamrcct(~vflagmissgene,:) ;
      caGeneNamenm = caGeneName(~vflagmissgene) ;

      %  Make PCA scatterplot
      %
      figure(1) ;
      clf ;
      titlecellstr = {{'PanCan Gene Expression, 2013' ...
                       [strct ', n = ' num2str(n) ' / ' num2str(nct)] ...
                        ['# No-missings genes = ' num2str(dnm) ' / ' num2str(d)]}} ;
      labelcellstr = {{'PC 1 Scores' 'PC 2 Scores' 'PC 3 Scores'}} ;
      savestr = ['PanCan1ip' num2str(ipart) 'PCA-NonMissGenes' strct] ;
      paramstruct = struct('npcadiradd',3, ...
                           'icolor',mcolorct, ...
                           'isubpopkde',1, ...
                           'ibelowdiag',0, ...
                           'titlecellstr',titlecellstr, ...
                           'labelcellstr',labelcellstr, ...
                           'savestr',savestr, ...
                           'iscreenwrite',1) ;
      scatplotSM(mdatamrcctnm,[],paramstruct) ;

      %  Make Mean MargDistPlot
      %
      figure(2) ;
      clf ;
      varnamecellstr = {caGeneNamenm} ;
      titlecellstr = {{'Summary Stats (means)' ...
                       [strct ', Non-Missing Genes'] ...
                       '7 smallest & 8 biggest'}} ;
      savestr = ['PanCan1ip' num2str(ipart) 'MargDistMean' strct] ;
      paramstruct = struct('istat',1, ...
                           'varnamecellstr',varnamecellstr, ...
                           'viplot',[(1:7), ((dnm - 7):dnm)], ...
                           'icolor',mcolorct, ...
                           'isubpopkde',1, ...
                           'datovlaymax',0.9, ...
                           'datovlaymin',0.4, ...
                           'titlecellstr',titlecellstr, ...
                           'savestr',savestr) ;
      MargDistPlotSM(mdatamrcctnm,paramstruct) ;


    end ;    %  of ict loop


  elseif ipart == 5 ;    %  Explore Subset COAD, READ, UCEC

    %  First do median Recentering of each gene
    %
    mdatamrc = [] ;
    for i = 1:d ;
      if (i / 1000) == floor(i / 1000) ;
         disp(['    Working on gene median ' num2str(i)]) ;
      end ;
      vdata = mdata(i,:) ;
      vflagmiss = isnan(vdata) ;
      med = median(vdata(~vflagmiss)) ;
      vdatamrc = vdata - med ;
      mdatamrc = [mdatamrc; vdatamrc] ;
    end ;
    disp('Finished Median Recentering') ;

    %  Find subsets
    %
    vflagss = vCOADflag | vREADflag | vUCECflag ;
    nss = sum(vflagss) ;
    disp(' ') ;
    disp(['  Check number of cases in subset is: ' num2str(nss)]) ;
    mdatamrcss = mdata(:,vflagss) ;
    mcolorss = mcolor(vflagss',:) ;
    vCOADflagss = vCOADflag(vflagss) ;
    vREADflagss = vREADflag(vflagss) ;
    vUCECflagss = vUCECflag(vflagss) ;

    %  Work only with genes with no missings
    %
    mflagmiss = isnan(mdatamrcss) ;
    vnmissgene = sum(mflagmiss,2) ;
    vflagmissgene = (vnmissgene > 0) ;
        %  one for genes with a missing
    dnm = sum(~vflagmissgene) ;
    disp(' ') ;
    disp(['  Check number of non-missing genes is: ' num2str(dnm)]) ;
    mdatanm = mdatamrcss(~vflagmissgene,:) ;
    caGeneNamenm = caGeneName(~vflagmissgene) ;

    %  Make PCA scatterplot
    %
    figure(1) ;
    clf ;
    titlecellstr = {{'PanCan Gene Expression, 2013' ...
                     'Types COAD, READ, UCEC only' ...
                     ['n = ' num2str(nss) ' / ' num2str(n)] ...
                     ['# No-missings genes = ' num2str(dnm) ' / ' num2str(d)]}} ;
    labelcellstr = {{'PC 1 Scores' 'PC 2 Scores' 'PC 3 Scores' 'PC 4 Scores'}} ;
    savestr = ['PanCan1ip' num2str(ipart) 'COADREADUCEC-PCA'] ;
    paramstruct = struct('npcadiradd',4, ...
                         'icolor',mcolorss, ...
                         'isubpopkde',1, ...
                         'ibelowdiag',0, ...
                         'titlecellstr',titlecellstr, ...
                         'labelcellstr',labelcellstr, ...
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    scatplotSM(mdatanm,[],paramstruct) ;

    %  Compute DWD directions
    %
    DWDUCECrestdirn = DWD2XQ(mdatanm(:,vUCECflagss),mdatanm(:,~vUCECflagss)) ;
    DWDCOADREADdirn = DWD2XQ(mdatanm(:,vCOADflagss),mdatanm(:,vREADflagss)) ;

    %  Make DWD scatterplot
    %
    figure(2) ;
    clf ;
    titlecellstr = {{'PanCan Gene Expression, 2013' ...
                     'Types COAD, READ, UCEC only' ...
                     ['n = ' num2str(nss) ' / ' num2str(n)] ...
                     ['# No-missings genes = ' num2str(dnm) ' / ' num2str(d)]}} ;
    labelcellstr = {{'DWD ECEC v Rest' 'DWD COAD v READ' 'OPC 1 Scores'}} ;
    savestr = ['PanCan1ip' num2str(ipart) 'COADREADUCEC-DWDs'] ;
    paramstruct = struct('npcadiradd',-1, ...
                         'icolor',mcolorss, ...
                         'isubpopkde',1, ...
                         'ibelowdiag',1, ...
                         'titlecellstr',titlecellstr, ...
                         'labelcellstr',labelcellstr, ...
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    scatplotSM(mdatanm,[DWDUCECrestdirn DWDCOADREADdirn],paramstruct) ;

    %  DiProPerm UCEC vs. rest
    %
    figure(3) ;
    clf ;
    icolor = [mcolorbase(11,:); mean([mcolorbase(2,:); mcolorbase(8,:)],1)] ;
    savestr = ['PanCan1ip' num2str(ipart) 'DiProPermUCECrest'] ;
    paramstruct = struct('nsim',100, ...
                         'nreport',20, ...
                         'icolor',icolor, ...
                         'title1str','DiProPerm UCEC v, Rest', ...
                         'title2str','UCEC-COAD-READ PanCan data', ...
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    DiProPermSM(mdatanm(:,vUCECflagss),mdatanm(:,~vUCECflagss),paramstruct) ;

    %  DiProPerm COAD vs. READ
    %
    figure(4) ;
    clf ;
    icolor = [mcolorbase(2,:); mcolorbase(8,:)] ;
    savestr = ['PanCan1ip' num2str(ipart) 'DiProPermCOADvREAD'] ;
    paramstruct = struct('nsim',100, ...
                         'nreport',20, ...
                         'icolor',icolor, ...
                         'title1str','DiProPerm COAD v, READ', ...
                         'title2str','UCEC-COAD-READ PanCan data', ...
                         'savestr',savestr, ...
                         'iscreenwrite',1) ;
    DiProPermSM(mdatanm(:,vCOADflagss),mdatanm(:,vREADflagss),paramstruct) ;

    %  Look at drivers of COAD vs. READ difference
    %
    titlestr = 'DWD Loadings, COAD vs. READ' ;
    savestr = ['PanCan1ip' num2str(ipart) 'DWDLoadingsCOADvREAD'] ;
    paramstruct = struct('isort',2, ...
                         'nshow',20, ...
                         'fontsize',15, ...
                         'titlestr',titlestr, ...
                         'ylabelstr','DWD Loadings', ...
                         'savestr',savestr) ;
    LabeledBarPlotSM(DWDCOADREADdirn,caGeneNamenm,paramstruct) ;


  elseif ipart == 6 ;    %  Classify each type versus the rest

    %  First do median Recentering of each gene
    %
    mdatamrc = [] ;
    for i = 1:d ;
      vdata = mdata(i,:) ;
      vflagmiss = isnan(vdata) ;
      med = median(vdata(~vflagmiss)) ;
      vdatamrc = vdata - med ;
      mdatamrc = [mdatamrc; vdatamrc] ;
    end ;
    disp('Finished Median Recentering') ;

    %  Kick out genes with missings
    %
    mflagmiss = isnan(mdatamrc) ;
    vnmissgene = sum(mflagmiss,2) ;
    vflagmissgene = (vnmissgene > 0) ;
        %  one for genes with a missing
    dnm = sum(~vflagmissgene) ;
    mdatamrcnm = mdatamrc(~vflagmissgene,:) ;
    caGeneNamenm = caGeneName(~vflagmissgene) ;

    for ict = 1:12 ;

      strct = caGeneNamebase{ict} ;
      vCTflag = mCTflag(ict,:) ;
      nct = sum(vCTflag) ;
%{
      mdatamrcct = mdatamrc(:,vCTflag) ;
%}
      vgray = [0.5 0.5 0.5] ;
      icolor = [mcolorbase(ict,:); vgray] ;

      mcolorct = ones(n,1) * vgray ;
      mcolorct(vCTflag,:) = mcolor(vCTflag,:) ;

      disp(' ') ;
      disp(['    Working on Cancer Type ' strct]) ;


      %  Make Marginal Overlap Plot    
      %
      figure(1) ;
      clf ;
      titlecellstr = {{['Marginal Overlap, T-pval'] ...
                       ['PanCan, 2013, ' strct ' vs. Rest'] ...
                       ['n = ' num2str(nct) ' / ' num2str(n) ...
                        'd = ' num2str(dnm) ' / ' num2str(d)]}} ;
      savestr = ['PanCan1ip' num2str(ipart) strct '-MargOverlapvRest'] ;
      paramstruct = struct('ioverlap',1, ...
                           'varnamecellstr',{caGeneNamenm}, ...
                           'icolor',icolor, ...
                           'isubpopkde',1, ...
                           'titlecellstr',titlecellstr, ...
                           'savestr',savestr) ;
      MargOverlapSM(mdatamrcnm,vCTflag,paramstruct) ;


      %  Compute DWD directions
      %
%      DWDdirn = DWD2XQ(mdatamrcnm(:,vCTflag),mdatamrcnm(:,~vCTflag)) ;
%{
      y = [ones(sum(vCTflag),1); -ones(sum(~vCTflag),1)] ; 
      [C,ddist] = penaltyParameter(mdatamrcnm,y,1); 
      DWDdirn = genDWDweighted([mdatamrcnm(:,vCTflag) mdatamrcnm(:,~vCTflag)],y,C,1) ;
          %  DWD Fast direction vector, pointing from 2nd group towards first
%}
      MDdirn = mean(mdatamrcnm(:,vCTflag),2) - mean(mdatamrcnm(:,~vCTflag),2) ;
      MDdirn = MDdirn / norm(MDdirn) ;

      %  Make DWD scatterplot
      %
      figure(2) ;
      clf ;
      titlecellstr = {{['PanCan, 2013, ' strct ' vs. Rest'] ...
                       ['n = ' num2str(nct) ' / ' num2str(n) ...
                        'd = ' num2str(dnm) ' / ' num2str(d)]}} ;
      labelcellstr = {{['MD '  strct ' v Rest'] 'OPC 1 Scores'}} ;
      savestr = ['PanCan1ip' num2str(ipart) strct '-MDvsRestnOPCA'] ;
      paramstruct = struct('npcadiradd',-1, ...
                           'icolor',mcolorct, ...
                           'isubpopkde',1, ...
                           'ibelowdiag',1, ...
                           'titlecellstr',titlecellstr, ...
                           'labelcellstr',labelcellstr, ...
                           'savestr',savestr, ...
                           'iscreenwrite',1) ;
%      scatplotSM(mdatamrcnm,DWDdirn,paramstruct) ;
      scatplotSM(mdatamrcnm,MDdirn,paramstruct) ;

      %  DiProPerm 
      %
      figure(3) ;
      clf ;
      icolor = [mcolorbase(ict,:); vgray] ;
      savestr = ['PanCan1ip' num2str(ipart) strct '-MDvRestDiProPerm'] ;
      paramstruct = struct('idir',2, ...
                           'nsim',100, ...
                           'nreport',5, ...
                           'icolor',icolor, ...
                           'title1str','PanCan 2013 data', ...
                           'title2str',['DiProPerm ' strct ' v, Rest'], ...
                           'savestr',savestr, ...
                           'iscreenwrite',1) ;
      DiProPermSMnew(mdatamrcnm(:,vCTflag),mdatamrcnm(:,~vCTflag),paramstruct) ;

      %  Look at drivers of difference
      %
      figure(4) ;
      clf ;
      titlestr = ['DWD Loadings, ' strct ' vs. Rest'] ;
      savestr = ['PanCan1ip' num2str(ipart) strct '-MDvRestLoadings'] ;
      paramstruct = struct('isort',2, ...
                           'nshow',20, ...
                           'fontsize',15, ...
                           'titlestr',titlestr, ...
                           'ylabelstr','MD Loadings', ...
                           'savestr',savestr) ;
      LabeledBarPlotSM(MDdirn,caGeneNamenm,paramstruct) ;


    end ;    %  of ict loop


  elseif ipart == 10 ;    %  Do better handling of missings, and store as .mat file

    mflagmiss = isnan(mdata) ;
    vnmissgene = sum(mflagmiss,2) ;
    vnmisscase = sum(mflagmiss,1) ;
    disp(' ') ;
    disp(['    # genes w/ missings = ' num2str(sum(vnmissgene > 0)) ...
                   '  out of ' num2str(d)]) ;
    disp(['    # cases w/ missings = ' num2str(sum(vnmisscase > 0)) ...
                   '  out of ' num2str(n)]) ;
    nmiss = num2str(sum(sum(mflagmiss,2))) ;
    disp(['       Total # missings = ' num2str(nmiss) ...
                   '  out of ' num2str(n * d)]) ;

    %  Map missings into 0s
    %
    mdatam0 = [] ;
    for j = 1:n ;
      vdata = mdata(:,j) ;
      vdata(mflagmiss(:,j)) = zeros(vnmisscase(j),1) ;
      mdatam0 = [mdatam0 vdata] ;
    end ;
    mdata = [] ;
    mflagmiss = [] ;
        %  to save space
    disp(' ') ;
    disp(['    Check All missings have been replaced:  ' ...
               num2str(sum(sum(isnan(mdatam0))))]) ;

    %  Check minimum values using MargDistPlot
    %
    figure(1) ;
    clf ;
    varnamecellstr = {caGeneName} ;
    titlecellstr = {{'Summary Stats (mins)' ...

                     'Raw RNAseq data' ...
                     'Missing Genes to 0'}} ;
    savestr = ['PanCan1ip' num2str(ipart) 'RawMargDistMin'] ;
    paramstruct = struct('istat',8, ...
                         'varnamecellstr',varnamecellstr, ...
                         'icolor',mcolor, ...
                         'isubpopkde',1, ...
                         'datovlaymax',0.9, ...
                         'datovlaymin',0.4, ...
                         'titlecellstr',titlecellstr, ...
                         'savestr',savestr) ;
    MargDistPlotSM(mdatam0,paramstruct) ;

    figure(2) ;
    clf ;
    varnamecellstr = {caGeneName} ;
    titlecellstr = {{'Summary Stats (mins)' ...
                     'Raw RNAseq data' ...
                     'Missing Genes to 0' ...
                     'Smallest 15'}} ;
    savestr = ['PanCan1ip' num2str(ipart) 'RawMargDistMinSmallest'] ;
    paramstruct = struct('istat',8, ...
                         'varnamecellstr',varnamecellstr, ...
                         'viplot',(1:15), ...
                         'icolor',mcolor, ...
                         'isubpopkde',1, ...
                         'datovlaymax',0.9, ...
                         'datovlaymin',0.4, ...
                         'titlecellstr',titlecellstr, ...
                         'savestr',savestr) ;
    MargDistPlotSM(mdatam0,paramstruct) ;
    
    figure(3) ;
    clf ;
    varnamecellstr = {caGeneName} ;
    titlecellstr = {{'Summary Stats (Num 0s)' ...
                     'Raw RNAseq data' ...
                     'Missing Genes to 0' ...
                     ['Largest, out of n = ' num2str(n)]}} ;
    savestr = ['PanCan1ip' num2str(ipart) 'RawMargDistNum0s'] ;
    paramstruct = struct('istat',13, ...
                         'varnamecellstr',varnamecellstr, ...
                         'viplot',((d - 14):d), ...
                         'icolor',mcolor, ...
                         'isubpopkde',1, ...
                         'datovlaymax',0.9, ...
                         'datovlaymin',0.4, ...
                         'titlecellstr',titlecellstr, ...
                         'savestr',savestr) ;
    MargDistPlotSM(mdatam0,paramstruct) ;

    figure(4) ;
    clf ;
    varnamecellstr = {caGeneName} ;
    titlecellstr = {{'Summary Stats (means)' ...
                     'Raw RNAseq data' ...
                     'Missing Genes to 0' ...
                     '7 smallest & 8 biggest'}} ;
    savestr = ['PanCan1ip' num2str(ipart) 'RawMargDistMeans'] ;
    paramstruct = struct('istat',1, ...
                         'varnamecellstr',varnamecellstr, ...
                         'viplot',[(1:7), ((d - 7):d)], ...
                         'icolor',mcolor, ...
                         'isubpopkde',1, ...
                         'datovlaymax',0.9, ...
                         'datovlaymin',0.4, ...
                         'titlecellstr',titlecellstr, ...
                         'savestr',savestr) ;
    MargDistPlotSM(mdatam0,paramstruct) ;


    %  Make Color Plot with Cancer Type Colors
    %
    figure(5) ;
    clf ;
    axis([0 1 0 1]) ;
    nct = 12 ;
    text(0.2,(nct + 1) / (nct + 2), ...
         'Type','Color','k','FontSize',24) ;
    text(0.4,(nct + 1) / (nct + 2), ...
         'n','Color','k','FontSize',24) ;
    text(0.6,(nct + 1) / (nct + 2), ...
         '#0s','Color','k','FontSize',24) ;
    text(0.8,(nct + 1) / (nct + 2), ...
         '#<0','Color','k','FontSize',24) ;
    for ict = 1:nct ;
      mflagE0CT = (mdatam0(:,mCTflag(ict,:)) == 0) ;
      mflagLT0CT = (mdatam0(:,mCTflag(ict,:)) < 0) ;
      text(0.2,(nct - ict + 1) / (nct + 2), ...
           caGeneNamebase{ict}, ...
           'Color',mcolorbase(ict,:), ...
           'FontSize',24) ;
      text(0.4,(nct - ict + 1) / (nct + 2), ...
           num2str(sum(mCTflag(ict,:))), ...
           'Color',mcolorbase(ict,:), ...
           'FontSize',24) ;
      text(0.6,(nct - ict + 1) / (nct + 2), ...
           num2str(sum(sum(mflagE0CT))), ...
           'Color',mcolorbase(ict,:), ...
           'FontSize',24) ;
      text(0.8,(nct - ict + 1) / (nct + 2), ...
           num2str(sum(sum(mflagLT0CT))), ...
           'Color',mcolorbase(ict,:), ...
           'FontSize',24) ;
    end ;
    orient landscape ;
    savestr = ['PanCan1ip' num2str(ipart) 'CancerTypeColors-NumMiss'] ;
    print('-dpsc2',savestr) ;


    %  Move negatives to 0 and do median recentering of each gene
    %
    mdatam0mrc = [] ;
    for i = 1:d ;
      if floor(i / 100) == (i / 100) ;
        disp(['        Median Recentering gene ' num2str(i)]) ;
      end ;
      vdatam0 = mdatam0(i,:) ;
      vdatam0 = max(vdatam0,0) ;
          %  replace negative values with 0
      vflag0 = (vdatam0 == 0) ;
      n0 = sum(vflag0) ;
      if n0 == n ;    %  Then all cases are 0
        med = 0 ;
      else ;
        med = median(vdatam0(~vflag0)) ;
      end ;
      vdatamrc = vdatam0 - med ;
      mdatam0mrc = [mdatam0mrc; vdatamrc] ;
    end ;
    disp('Finished Median Recentering') ;
    mdatam0 = [] ;
        %  to save space

    %  Check results with MargDistPlot, mean
    %
    figure(6) ;
    clf ;
    varnamecellstr = {caGeneName} ;
    titlecellstr = {{'Summary Stats (means)' ...
                     'Normalized data' ...
                     '7 smallest & 8 biggest'}} ;
    savestr = ['PanCan1ip' num2str(ipart) 'NormalizedMargDistMeans'] ;
    paramstruct = struct('istat',1, ...
                         'varnamecellstr',varnamecellstr, ...
                         'viplot',[(1:7), ((d - 7):d)], ...
                         'icolor',mcolor, ...
                         'isubpopkde',1, ...
                         'datovlaymax',0.9, ...
                         'datovlaymin',0.4, ...
                         'titlecellstr',titlecellstr, ...
                         'savestr',savestr) ;
    MargDistPlotSM(mdatam0mrc,paramstruct) ;

    mdatan = mdatam0mrc ;
    mdatam0mrc = [] ;
        %  to save space

    disp(' ') ;
    disp('    Start .mat file save') ;
    save(savefileAstr,'mdatan') ;
    disp('    Finished .mat file save') ;
    disp(' ') ;

  else ;    %  Read in normalized version of data
            %  Missings and negatives properly handled
            %      and median recentered

    disp(' ') ;
    disp('    Start .mat file read') ;
    load(savefileAstr) ;
    %    Loads Variable:
    %        mdatan:         Matrix of Gene Expression

    %    Note most other variables loaded above
    mdata = [] ;
        %  To save space, this was loaded above, but not needed


    if ipart == 11 ;    %  PCA Scatterplots

      %  PCA scatterplot, all together
      %
      figure(2) ;
      clf ;
      titlecellstr = {{'PanCan Gene Expression, 2013' ...
                       [num2str(size(mcolorbase,1)) 'Cancer Types'] ...
                       ['n = ' num2str(n)] ...
                       ['# genes, d = ' num2str(d)]}} ;
      labelcellstr = {{'PC 1 Scores' 'PC 2 Scores' 'PC 3 Scores' 'PC 4 Scores'}} ;
      savestr = ['PanCan1ip' num2str(ipart) 'AllDataPCA'] ;
      paramstruct = struct('npcadiradd',4, ...
                           'icolor',mcolor, ...
                           'isubpopkde',1, ...
                           'ibelowdiag',0, ...
                           'titlecellstr',titlecellstr, ...
                           'labelcellstr',labelcellstr, ...
                           'savestr',savestr, ...
                           'iscreenwrite',1) ;
      scatplotSM(mdatan,[],paramstruct) ;

      %  Make some zoomed in plots
      %
      paramstruct = struct('npc',4, ...
                           'iorient',3, ...
                           'iscreenwrite',1, ...
                           'viout',[0 0 0 0 1]) ;
      outstruct = pcaSM(mdatan,paramstruct) ;
      mpc = getfield(outstruct,'mpc') ;
          %  4 x n matrix of PC scores

      figure(3) ;
      clf ;
      savestr = ['PanCan1ip' num2str(ipart) 'AllData-PC1v2'] ;
      paramstruct = struct('icolor',mcolor, ...
                           'titlestr','Zoomed in PanCan All', ...
                           'xlabelstr','PC 2 Scores', ...
                           'ylabelstr','PC 1 Scores', ...
                           'ifigure',3, ...
                           'savestr',savestr, ...
                           'iscreenwrite',1) ;
      projplot2SM(mpc,[[0; 1; 0; 0],[1; 0; 0; 0]],paramstruct) ;

      figure(4) ;
      clf ;
      savestr = ['PanCan1ip' num2str(ipart) 'AllData-PC2v3'] ;
      paramstruct = struct('icolor',mcolor, ...
                           'titlestr','Zoomed in PanCan All', ...
                           'xlabelstr','PC 3 Scores', ...
                           'ylabelstr','PC 2 Scores', ...
                           'ifigure',4, ...
                           'savestr',savestr, ...
                           'iscreenwrite',1) ;
      projplot2SM(mpc,[[0; 0; 1; 0],[0; 1; 0; 0]],paramstruct) ;

      figure(5) ;
      clf ;
      savestr = ['PanCan1ip' num2str(ipart) 'AllData-PC3v4'] ;
      paramstruct = struct('icolor',mcolor, ...
                           'titlestr','Zoomed in PanCan All', ...
                           'xlabelstr','PC 4 Scores', ...
                           'ylabelstr','PC 3 Scores', ...
                           'ifigure',5, ...
                           'savestr',savestr, ...
                           'iscreenwrite',1) ;
      projplot2SM(mpc,[[0; 0; 0; 1],[0; 0; 1; 0]],paramstruct) ;


    elseif ipart == 12 ;    %  Individual Mean MargDist plots & PCAs

      for ict = 1:12 ;

        strct = caGeneNamebase{ict} ;
        vflagCT = mCTflag(ict,:) ;
        mdatanct = mdatan(:,vflagCT) ;
        mcolorct = mcolor(vflagCT,:) ;
        nct = sum(vflagCT) ;

        disp(['    Working on Cancer Type ' strct]) ;

        %  Make PCA scatterplot
        %
        figure(1) ;
        clf ;
        titlecellstr = {{'PanCan Gene Expression, 2013' ...
                         [strct ', n = ' num2str(nct) ' / ' num2str(n)] ...
                         ['Normalized Data, d = ' num2str(d)]}} ;
        labelcellstr = {{'PC 1 Scores' 'PC 2 Scores' 'PC 3 Scores'}} ;
        savestr = ['PanCan1ip' num2str(ipart) strct 'PCA'] ;
        paramstruct = struct('npcadiradd',3, ...
                             'icolor',mcolorct, ...
                             'ibelowdiag',0, ...
                             'titlecellstr',titlecellstr, ...
                             'labelcellstr',labelcellstr, ...
                             'savestr',savestr, ...
                             'iscreenwrite',1) ;
        scatplotSM(mdatanct,[],paramstruct) ;

        %  Make Mean MargDistPlot
        %
        figure(2) ;
        clf ;
        varnamecellstr = {caGeneName} ;
        titlecellstr = {{'Summary Stats (means)' ...
                         strct ...
                         'Full data median centered' ...
                         '7 smallest & 8 biggest'}} ;
        savestr = ['PanCan1ip' num2str(ipart) strct 'MargDistMean'] ;
        paramstruct = struct('istat',1, ...
                             'varnamecellstr',varnamecellstr, ...
                             'viplot',[(1:7), ((d - 7):d)], ...
                             'icolor',mcolorct, ...
                             'datovlaymax',0.9, ...
                             'datovlaymin',0.4, ...
                             'titlecellstr',titlecellstr, ...
                             'savestr',savestr) ;
        MargDistPlotSM(mdatanct,paramstruct) ;


      end ;    %  of ict loop


    elseif ipart == 13 ;    %  Explore Subset COAD, READ, UCEC

%{
      %  First do median Recentering of each gene
      %
      mdatamrc = [] ;
      for i = 1:d ;
        if (i / 1000) == floor(i / 1000) ;
           disp(['    Working on gene median ' num2str(i)]) ;
        end ;
        vdata = mdata(i,:) ;
        vflagmiss = isnan(vdata) ;
        med = median(vdata(~vflagmiss)) ;
        vdatamrc = vdata - med ;
        mdatamrc = [mdatamrc; vdatamrc] ;
      end ;
      disp('Finished Median Recentering') ;
%}

      %  Find subsets
      %
      vflagss = vCOADflag | vREADflag | vUCECflag ;
      nss = sum(vflagss) ;
      disp(' ') ;
      disp(['  Check number of cases in subset is: ' num2str(nss)]) ;
      mdatanss = mdatan(:,vflagss) ;
      mcolorss = mcolor(vflagss',:) ;
      vCOADflagss = vCOADflag(vflagss) ;
      vREADflagss = vREADflag(vflagss) ;
      vUCECflagss = vUCECflag(vflagss) ;

      %  Make PCA scatterplot
      %
      figure(1) ;
      clf ;
      titlecellstr = {{'PanCan Gene Expression, 2013' ...
                       'Types COAD, READ, UCEC only' ...
                       ['n = ' num2str(nss) ' / ' num2str(n) ', Normalized'] ...
                       ['# genes, d = ' num2str(d)]}} ;
      labelcellstr = {{'PC 1 Scores' 'PC 2 Scores' 'PC 3 Scores' 'PC 4 Scores'}} ;
      savestr = ['PanCan1ip' num2str(ipart) 'COADREADUCEC-PCA'] ;
      paramstruct = struct('npcadiradd',4, ...
                           'icolor',mcolorss, ...
                           'isubpopkde',1, ...
                           'ibelowdiag',0, ...
                           'titlecellstr',titlecellstr, ...
                           'labelcellstr',labelcellstr, ...
                           'savestr',savestr, ...
                           'iscreenwrite',1) ;
      scatplotSM(mdatanss,[],paramstruct) ;

      %  Compute DWD directions
      %
      DWDUCECrestdirn = DWD2XQ(mdatanss(:,vUCECflagss),mdatanss(:,~vUCECflagss)) ;
      DWDCOADREADdirn = DWD2XQ(mdatanss(:,vCOADflagss),mdatanss(:,vREADflagss)) ;

      %  Make DWD scatterplot
      %
      figure(2) ;
      clf ;
      titlecellstr = {{'PanCan Gene Expression, 2013' ...
                       'Types COAD, READ, UCEC only' ...
                       ['n = ' num2str(nss) ' / ' num2str(n) ', Normalized'] ...
                       ['# genes, d = ' num2str(d)]}} ;
      labelcellstr = {{'DWD ECEC v Rest' 'DWD COAD v READ' 'OPC 1 Scores'}} ;
      savestr = ['PanCan1ip' num2str(ipart) 'COADREADUCEC-DWDs'] ;
      paramstruct = struct('npcadiradd',-1, ...
                           'icolor',mcolorss, ...
                           'isubpopkde',1, ...
                           'ibelowdiag',1, ...
                           'titlecellstr',titlecellstr, ...
                           'labelcellstr',labelcellstr, ...
                           'savestr',savestr, ...
                           'iscreenwrite',1) ;
      scatplotSM(mdatanss,[DWDUCECrestdirn DWDCOADREADdirn],paramstruct) ;

      %  DiProPerm UCEC vs. rest
      %
      figure(3) ;
      clf ;
      icolor = [mcolorbase(11,:); mean([mcolorbase(2,:); mcolorbase(8,:)],1)] ;
      savestr = ['PanCan1ip' num2str(ipart) 'DiProPermUCECrest'] ;
      paramstruct = struct('nsim',100, ...
                           'nreport',20, ...
                           'icolor',icolor, ...
                           'title1str','DiProPerm UCEC v, Rest', ...
                           'title2str','UCEC-COAD-READ PanCan data', ...
                           'savestr',savestr, ...
                           'iscreenwrite',1) ;
      DiProPermSM(mdatanss(:,vUCECflagss),mdatanss(:,~vUCECflagss),paramstruct) ;

      %  DiProPerm COAD vs. READ
      %
      figure(4) ;
      clf ;
      icolor = [mcolorbase(2,:); mcolorbase(8,:)] ;
      savestr = ['PanCan1ip' num2str(ipart) 'DiProPermCOADvREAD'] ;
      paramstruct = struct('nsim',100, ...
                           'nreport',20, ...
                           'icolor',icolor, ...
                           'title1str','DiProPerm COAD v, READ', ...
                           'title2str','UCEC-COAD-READ PanCan data', ...
                           'savestr',savestr, ...
                           'iscreenwrite',1) ;
      DiProPermSM(mdatanss(:,vCOADflagss),mdatanss(:,vREADflagss),paramstruct) ;

      %  Look at drivers of COAD vs. READ difference
      %
      figure(5) ;
      clf ;
      titlestr = 'DWD Loadings, COAD vs. READ' ;
      savestr = ['PanCan1ip' num2str(ipart) 'DWDLoadingsCOADvREAD'] ;
      paramstruct = struct('isort',2, ...
                           'nshow',20, ...
                           'fontsize',15, ...
                           'titlestr',titlestr, ...
                           'ylabelstr','DWD Loadings', ...
                           'savestr',savestr) ;
      LabeledBarPlotSM(DWDCOADREADdirn,caGeneName,paramstruct) ;


    elseif ipart == 14 ;    %  Classify each type versus the rest

      for ict = 1:12 ;

        strct = caGeneNamebase{ict} ;
        vCTflag = mCTflag(ict,:) ;
        nct = sum(vCTflag) ;

        vgray = [0.5 0.5 0.5] ;
        icolor = [mcolorbase(ict,:); vgray] ;

        mcolorct = ones(n,1) * vgray ;
        mcolorct(vCTflag,:) = mcolor(vCTflag,:) ;

        disp(' ') ;
        disp(['    Working on Cancer Type ' strct]) ;

        %  Make Marginal Overlap Plot    
        %
        figure(1) ;
        clf ;
        titlecellstr = {{['Marginal Overlap, T-pval'] ...
                         ['PanCan, 2013, ' strct ' vs. Rest'] ...
                         ['n = ' num2str(nct) ' / ' num2str(n) ...
                          'd = ' num2str(d)]}} ;
        savestr = ['PanCan1ip' num2str(ipart) strct '-MargOverlapvRest'] ;
        paramstruct = struct('ioverlap',1, ...
                             'varnamecellstr',{caGeneName}, ...
                             'icolor',icolor, ...
                             'isubpopkde',1, ...
                             'titlecellstr',titlecellstr, ...
                             'savestr',savestr) ;
        MargOverlapSM(mdatan,vCTflag,paramstruct) ;


        %  Compute DWD directions
        %
  %      DWDdirn = DWD2XQ(mdatamrcnm(:,vCTflag),mdatamrcnm(:,~vCTflag)) ;
  %{
        y = [ones(sum(vCTflag),1); -ones(sum(~vCTflag),1)] ; 
        [C,ddist] = penaltyParameter(mdatamrcnm,y,1); 
        DWDdirn = genDWDweighted([mdatamrcnm(:,vCTflag) mdatamrcnm(:,~vCTflag)],y,C,1) ;
            %  DWD Fast direction vector, pointing from 2nd group towards first
  %}
        MDdirn = mean(mdatan(:,vCTflag),2) - mean(mdatan(:,~vCTflag),2) ;
        MDdirn = MDdirn / norm(MDdirn) ;

        %  Make MD scatterplot
        %
        figure(2) ;
        clf ;
        titlecellstr = {{['PanCan, 2013, ' strct ' vs. Rest'] ...
                         ['n = ' num2str(nct) ' / ' num2str(n) ...
                          'd = ' num2str(d)]}} ;
        labelcellstr = {{['MD '  strct ' v Rest'] 'OPC 1 Scores'}} ;
        savestr = ['PanCan1ip' num2str(ipart) strct '-MDvsRestnOPCA'] ;
        paramstruct = struct('npcadiradd',-1, ...
                             'icolor',mcolorct, ...
                             'isubpopkde',1, ...
                             'ibelowdiag',1, ...
                             'titlecellstr',titlecellstr, ...
                             'labelcellstr',labelcellstr, ...
                             'savestr',savestr, ...
                             'iscreenwrite',1) ;
  %      scatplotSM(mdatamrcnm,DWDdirn,paramstruct) ;
        scatplotSM(mdatan,MDdirn,paramstruct) ;

        %  DiProPerm 
        %
        figure(3) ;
        clf ;
        icolor = [mcolorbase(ict,:); vgray] ;
        savestr = ['PanCan1ip' num2str(ipart) strct '-MDvRestDiProPerm'] ;
        paramstruct = struct('idir',2, ...
                             'nsim',100, ...
                             'nreport',5, ...
                             'icolor',icolor, ...
                             'title1str','PanCan 2013 data', ...
                             'title2str',['DiProPerm ' strct ' v, Rest'], ...
                             'savestr',savestr, ...
                             'iscreenwrite',1) ;
        DiProPermSMnew(mdatan(:,vCTflag),mdatan(:,~vCTflag),paramstruct) ;

        %  Look at drivers of difference
        %
        figure(4) ;
        clf ;
        titlestr = ['MD Loadings, ' strct ' vs. Rest'] ;
        savestr = ['PanCan1ip' num2str(ipart) strct '-MDvRestLoadings'] ;
        paramstruct = struct('isort',2, ...
                             'nshow',20, ...
                             'fontsize',15, ...
                             'titlestr',titlestr, ...
                             'ylabelstr','MD Loadings', ...
                             'savestr',savestr) ;
        LabeledBarPlotSM(MDdirn,caGeneName,paramstruct) ;


      end ;    %  of ict loop


    elseif ipart == 15 ;    %  MargOverlap, AUC-ROC, each versus the rest

      for ict = 1:12 ;

        strct = caGeneNamebase{ict} ;
        vCTflag = mCTflag(ict,:) ;
        nct = sum(vCTflag) ;

        vgray = [0.5 0.5 0.5] ;
        icolor = [mcolorbase(ict,:); vgray] ;

        mcolorct = ones(n,1) * vgray ;
        mcolorct(vCTflag,:) = mcolor(vCTflag,:) ;

        disp(' ') ;
        disp(['    Working on Cancer Type ' strct]) ;

        %  Make Marginal Overlap Plot    
        %
        figure(1) ;
        clf ;
        titlecellstr = {{['Marginal Overlap, AUC-ROC'] ...
                         ['PanCan, 2013, ' strct ' vs. Rest'] ...
                         ['n = ' num2str(nct) ' / ' num2str(n) ...
                          'd = ' num2str(d)]}} ;
        savestr = ['PanCan1ip' num2str(ipart) strct '-AUCMargOverlapvRest'] ;
        paramstruct = struct('ioverlap',5, ...
                             'varnamecellstr',{caGeneName}, ...
                             'icolor',icolor, ...
                             'isubpopkde',1, ...
                             'titlecellstr',titlecellstr, ...
                             'savestr',savestr) ;
        MargOverlapSM(mdatan,vCTflag,paramstruct) ;


      end ;    %  of ict loop


    elseif ipart == 16 ;    %  Explore LUAD vs. LUSC





    elseif ipart == 20 ;    %  Pairwise classifications


    end ;    %  of inner-inner ipart if-block


  end ;    %  of inner ipart if-block


end ;    %  of outer ipart if-block



