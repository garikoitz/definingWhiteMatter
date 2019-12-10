% add these to your path:
%addpath(genpath('/data/localhome/lagudiez/scitran'))
%addpath(genpath('/data/localhome/lagudiez/JSONio'))
%addpath(genpath('/data/localhome/lagudiez/definingWhiteMatter'))
%addpath(genpath('/data/localhome/lagudiez/vistasoft'))
%addpath(genpath('/data/localhome/lagudiez/AFQ'))



%% 
MAINDIR   = '/share/wandell/users/glerma/TESTDATA/DefiningWMtractography/lucas/';
doDetProb = {'det','prob'};  % 'det', 'prob' or both
fmm       = '100';  % Tractogram filtered to fibers shorter than, in mm
useCortex = true;   % If we want to use cortex as ROI

% Dictionary for ROI 2 TRACT correspondence
R2T = containers.Map();
R2T(   'L_ATR') = {'ATR_roi1_L.nii.gz','ATR_roi2_L.nii.gz' };
R2T(   'R_ATR') = {'ATR_roi1_R.nii.gz','ATR_roi2_R.nii.gz' };
R2T(   'L_CGC') = {'CGC_roi1_L.nii.gz','CGC_roi2_L.nii.gz' };
R2T(   'R_CGC') = {'CGC_roi1_R.nii.gz','CGC_roi2_R.nii.gz' };
R2T(   'L_CST') = {'CST_roi1_L.nii.gz','CST_roi2_L.nii.gz' };
R2T(   'R_CST') = {'CST_roi1_R.nii.gz','CST_roi2_R.nii.gz' };
R2T(   'FA'   ) = {'FA_L.nii.gz'      ,'FA_R.nii.gz'       };
R2T(   'FP'   ) = {'FP_L.nii.gz'      ,'FP_R.nii.gz'       };
R2T(   'L_HCC') = {'HCC_roi1_L.nii.gz','HCC_roi2_L.nii.gz' };
R2T(   'R_HCC') = {'HCC_roi1_R.nii.gz','HCC_roi2_R.nii.gz' };
R2T(   'L_IFO') = {'IFO_roi1_L.nii.gz','IFO_roi2_L.nii.gz' };
R2T(   'R_IFO') = {'IFO_roi1_R.nii.gz','IFO_roi2_R.nii.gz' };
R2T(   'L_ILF') = {'ILF_roi1_L.nii.gz','ILF_roi2_L.nii.gz' };
R2T(   'R_ILF') = {'ILF_roi1_R.nii.gz','ILF_roi2_R.nii.gz' };
R2T(   'L_SLF') = {'SLF_roi1_L.nii.gz','SLF_roi2_L.nii.gz' };
R2T(   'R_SLF') = {'SLF_roi1_R.nii.gz','SLF_roi2_R.nii.gz' };
R2T(   'L_AF' ) = {'SLF_roi1_L.nii.gz','SLFt_roi2_L.nii.gz'};
R2T(   'R_AF' ) = {'SLF_roi1_R.nii.gz','SLFt_roi2_R.nii.gz'};
R2T(   'L_UNC') = {'UNC_roi1_L.nii.gz','UNC_roi2_L.nii.gz' };
R2T(   'R_UNC') = {'UNC_roi1_R.nii.gz','UNC_roi2_R.nii.gz' };

% Read all the subjects
dirlist = dir(fullfile(MAINDIR, 'ROIs/s*'));
subs = {dirlist.name};
disp(subs);

structural = fullfile(MAINDIR, "Structural");
tractsoutdir = fullfile(MAINDIR, "tracts");
ROIdir = fullfile(MAINDIR, "ROIs");
makeSureDir(tractsoutdir);

for i=1 : length(subs)
   s = subs{i}; 
   disp(s);
   makeSureDir(fullfile(tractsoutdir, s));
   makeSureDir(fullfile(tractsoutdir, s, strcat('tracts-',fmm)))
   makeSureDir(fullfile(tractsoutdir, s, strcat('tracts-',fmm),'det'))
   makeSureDir(fullfile(tractsoutdir, s, strcat('tracts-',fmm),'prob'))
end    

tcks = fullfile(MAINDIR, join(["TCK-Streamlines-",fmm,"mm"], ''));

% Do for all subjects
for i=1 : length(subs)
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
   cmd          = join(['maskfilter -npass 3 ',binaparc,' dilate ',binaparcdil],'');
   if ~isfile(binaparcdil)
       spres = system(cmd);
   end
   
   % Select what we want to use
   if useCortex
       Cortex = join([' -include ',binaparc,' '],''); % Can be ' '
   else
       Cortex = ' ';
   end
   
   % Get path to the det and prob tractograms
   tractogram = containers.Map();
   tractogram('det')  = fullfile(tcks, s, join(['tracking-deterministic-',fmm,'mm.tck'], ''));
   tractogram('prob') = fullfile(tcks, s, join(['tracking-probabilistic-',fmm,'mm.tck'], ''));
   
   % Obtain tracts dict and loop them
   tracts = keys(R2T);
   disp(tractogram('prob'))    
   for j=1 : length(tracts)
       tr = tracts{j};
       % Read required ROIs
       v = R2T(tr);
       ROI1   = fullfile(ROIdir, s,'ROIs', v{1});
       ROI2   = fullfile(ROIdir, s,'ROIs', v{2});
       
       % Run the mrtrix code
       % Set options
       exclude   = ' '; % can be ' ' or '-exclude roi.nii.gz '
       maxlength = ' '; % can be ' ' or '-maxlength xx '
       minlength = ' '; % can be ' ' or '-minlength xx ' 
       force     = ' -force '; % can be ' ' or '-force '
       
       % det
       for k=1 : length(doDetProb)
           dodp = doDetProb{k};
           tracks_in  = tractogram(dodp);
           
           if useCortex
                tracks_out = fullfile(tractsoutdir, s, strcat('tracts-',fmm), dodp, strcat(tr,'_cortex.tck'));
           else
                tracks_out = fullfile(tractsoutdir, s, strcat('tracts-',fmm), dodp, strcat(tr,'.tck'));
           end
           disp(tracks_out)
           
           cmd = join(['tckedit -include ' , ROI1, ' -include ', ROI2, ' ', ...
                             Cortex, exclude, maxlength, minlength, force, ...
                             tracks_in, ' ', tracks_out],'');
           spres = system(cmd);
       end
   end
end



disp('FINISHED waypoints.m')

%%


function makeSureDir(dirpath)
    % shorter alternative: if ~exist(dirpath,'dir') mkdir(dirpath); end
    if dirpath(end) ~= '/', dirpath = dirpath+'/'; end
    if (exist(dirpath, 'dir') == 0), mkdir(dirpath); end
end