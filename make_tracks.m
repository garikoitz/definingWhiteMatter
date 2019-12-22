% add these to your path:
%addpath(genpath('/data/localhome/lagudiez/definingWhiteMatter'))
%addpath(genpath('/data/localhome/lagudiez/vistasoft'))
%addpath(genpath('/data/localhome/lagudiez/AFQ'))
%addpath(genpath('/data/localhome/lagudiez/freesurfer_mrtrix_afni_matlab_tools'))


addpath(genpath('/Users/glerma/soft/definingWhiteMatter'))
MAINDIR   = '/Users/glerma/Documents/STANFORD_PROJECTS/DefiningWM';

% Dictionary for ROI 2 TRACT correspondence
R2T = containers.Map();
R2T('L_AF' ) = {'SLF_roi1_L.nii.gz', 'SLFt_roi2_L.nii.gz','', 'prob', '150', 3, 3, 1, true };
R2T('R_AF' ) = {'SLF_roi1_R.nii.gz', 'SLFt_roi2_R.nii.gz','', 'prob', '150', 3, 3, 1, true };
R2T('L_SLF') = {'SLF_roi1_L.nii.gz', 'SLF_roi2_L.nii.gz', '', 'prob', '100', 3, 3, 1, true };
R2T('R_SLF') = {'SLF_roi1_R.nii.gz', 'SLF_roi2_R.nii.gz', '', 'prob', '100', 3, 3, 1, true };
R2T('L_CST') = {'CST_roi1_L.nii.gz', 'CST_roi2_L.nii.gz', '', 'prob', '150', 3, 3, 1, true };
R2T('R_CST') = {'CST_roi1_R.nii.gz', 'CST_roi2_R.nii.gz', '', 'prob', '150', 3, 3, 1, true };
R2T('L_CGC') = {'CGC_roi1_L.nii.gz', 'CGC_roi2_L.nii.gz', '', 'prob', '150', 3, 3, 1, true };
R2T('R_CGC') = {'CGC_roi1_R.nii.gz', 'CGC_roi2_R.nii.gz', '', 'prob', '150', 3, 3, 1, true }
R2T('L_UNC') = {'UNC_roi1_L.nii.gz', 'UNC_roi2_L.nii.gz', '', 'prob', '100', 3, 3, 4, false }; %
R2T('R_UNC') = {'UNC_roi1_R.nii.gz', 'UNC_roi2_R.nii.gz', '', 'prob', '100', 3, 3, 4, false }; %
R2T('L_IFO') = {'IFO_roi1_L.nii.gz', 'IFO_roi2_L.nii.gz', 'IFO_roi3_L.nii.gz', 'prob', '200', 3, 3, 1, true};
R2T('R_IFO') = {'IFO_roi1_R.nii.gz', 'IFO_roi2_R.nii.gz', 'IFO_roi3_R.nii.gz', 'prob', '200', 3, 3, 1, true};
R2T('L_ILF') = {'ILF_roi1_L.nii.gz', 'ILF_roi2_L.nii.gz', 'ILF_roi3_L.nii.gz', 'prob', '150', 3, 3, 2, true }; %
R2T('R_ILF') = {'ILF_roi1_R.nii.gz', 'ILF_roi2_R.nii.gz', 'ILF_roi3_R.nii.gz', 'prob', '150', 3, 3, 2, true }; %


% Read all the subjects


dirlist = dir(fullfile(MAINDIR, 'ROIs/s*'));
subs = {dirlist.name};

structural = fullfile(MAINDIR, "Structural");
tractsoutdir = fullfile(MAINDIR, "tracts");
ROIdir = fullfile(MAINDIR, "ROIs");

if ~exist(MAINDIR,'dir') mkdir(MAINDIR); end
if ~exist(structural,'dir') mkdir(structural); end
if ~exist(tractsoutdir,'dir') mkdir(tractsoutdir); end
if ~exist(ROIdir,'dir') mkdir(ROIdir); end
for i=1 : length(subs)
   s = subs{i}; 
   if ~exist(fullfile(tractsoutdir, s),'dir') mkdir(fullfile(tractsoutdir, s)); end

end    

tcks = containers.Map();
tcks('100') = fullfile(MAINDIR, join(["TCK-Streamlines-",'100',"mm"],''));
tcks('150') = fullfile(MAINDIR, join(["TCK-Streamlines-",'150',"mm"],''));
tcks('200') = fullfile(MAINDIR, join(["TCK-Streamlines-",'200',"mm"],''));


% Do for all subjects
for i=1: length(subs)
   s = subs{i}; 
   % See if aparcaseg is binarized
   aparc    = fullfile(structural, s,'T1w','aparc+aseg.nii.gz ');
   binaparc = fullfile(structural, s,'T1w', 'segmentation.nii.gz');
   cmd      = join(['mrthreshold  -quiet -abs 1000 ',aparc,binaparc],' ');
   if ~isfile(binaparc)
       spres = system(cmd);
   end
   
   % And now dilate it so that we are sure that the fibers arrive
   binaparcdil  = fullfile(structural, s, 'T1w', 'segmentation_dilated.nii.gz');
   cmd          = join(['maskfilter -npass 3 ', binaparc,' dilate ', binaparcdil],'');
   if ~isfile(binaparcdil)
       spres = system(cmd);
   end
   
  
   % Get path to the det and prob tractograms
   tractogram = containers.Map();
   tractogram('det100')  = fullfile(tcks('100'), s, join(['tracking-deterministic-','100','mm.tck']));
   tractogram('prob100') = fullfile(tcks('100'), s, join(['tracking-probabilistic-','100','mm.tck']));
   tractogram('det150')  = fullfile(tcks('150'), s, join(['tracking-deterministic-','150','mm.tck']));
   tractogram('prob150') = fullfile(tcks('150'), s, join(['tracking-probabilistic-','150','mm.tck']));
   tractogram('det200')  = fullfile(tcks('200'), s, join(['tracking-deterministic-','200','mm.tck']));
   tractogram('prob200') = fullfile(tcks('200'), s, join(['tracking-probabilistic-','200','mm.tck']));
   
   
   % Obtain tracts dict and loop them
   tracts = keys(R2T);
   
   for j=1 : length(tracts)
       tr = tracts{j};
       % Read required ROIs
       v = R2T(tr);
       include_ROI1   = join([' -include ', fullfile(ROIdir, s,'ROIs', v{1})]); % include ROI
       include_ROI2   = join([' -include ', fullfile(ROIdir, s,'ROIs', v{2})]); % include ROI
       include_ROI3   = join([' -include ', fullfile(ROIdir, s,'ROIs', v{3})]); % include ROI
       if isempty(v{3})
          include_ROI3 = ' '; 
       end
       detprob= v{4};
       fmm    = v{5};
       maxDist= v{6};
       maxLen = v{7};
       maxIter= v{8};
       useCortex=v{9};
       % Select what we want to use
       if useCortex
           Cortex = join([' -include ', binaparc, ' '],''); % Can be ' '
       else
           Cortex = ' ';
       end
       
       % Run the mrtrix code
       % Set options
       include   = join([include_ROI1, include_ROI2, include_ROI3]);
       exclude   = ' '; % can be ' ' or '-exclude roi.nii.gz '
       maxlength = ' '; % can be ' ' or '-maxlength xx ' >>> tractogram already selected
       minlength = ' '; % can be ' ' or '-minlength xx ' 
       force     = ' -force '; % can be ' ' or '-force '
       
       tracks_in  = tractogram(join([detprob,fmm]));

       if useCortex
            tracks_out = fullfile(tractsoutdir, s, strcat(tr,'_cortex.tck'));
       else
            tracks_out = fullfile(tractsoutdir, s, strcat(tr,'.tck'));
       end
       
       cmd = join(['tckedit -quiet ', exclude, ' ', include, Cortex, maxlength, minlength, force, ...
                         tracks_in, ' ', tracks_out],'');
  
       
       % disp(cmd)
       spres = system(cmd);
       
       fileattrib(tracks_in, '+w +x') % make it readable and writeable
       fileattrib(tracks_out, '+w +x') % make it readable and writeable
       
       
       %%
       if useCortex
            tracks_original = fullfile(tractsoutdir, s, strcat(tr,'_cortex.tck'));
            tracks_cleaned = fullfile(tractsoutdir, s, join(['cleaned_',tr,'_cortex.tck'],''));
            tracks_image = fullfile(tractsoutdir, s, join(['cleaned_',tr,'_cortex.png'],''));
            tracks_image_original = fullfile(tractsoutdir, s, join([tr,'_cortex.png'],''));

       else
            tracks_original = fullfile(tractsoutdir, s, strcat(tr,'.tck'));
            tracks_cleaned = fullfile(tractsoutdir, s, join(['cleaned_',tr,'.tck'],''));
            tracks_image = fullfile(tractsoutdir, s, join(['cleaned_',tr,'.png'],''));
            tracks_image_original = fullfile(tractsoutdir, s, join([tr,'.png'],''));
       end


       track = fgRead(tracks_original);

       if size(track.fibers{1},2) < 1
           disp(join(['There are no fibers for track ', tr, ' Skipping.'],''));
           continue
       end


       clean_track = AFQ_removeFiberOutliers(track,maxDist,maxLen,100,'median',1,maxIter);

       fgWrite(clean_track, tracks_cleaned,'tck');
       fileattrib(tracks_cleaned, '+w +x') % make it readable and writeable


   end
end

disp('Finished')