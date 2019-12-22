# -*- coding: utf-8 -*-
import os
import glob as glob
import subprocess as sp
join      = os.path.join

# Stablish directories and options
SUBJDIR   = '/share/wandell/users/glerma/TESTDATA/DefiningWMtractography/lucas/'
doDetProb = ['det','prob']  # 'det', 'prob' or both
fmm       = '100'  # Tractogram filtered to fibers shorter than, in mm
useCortex = True   # If we want to use cortex as ROI

# Dictionary for ROI 2 TRACT correspondence
R2T = {'L_ATR':('ATR_roi1_L.nii.gz','ATR_roi2_L.nii.gz' ),
       'R_ATR':('ATR_roi1_R.nii.gz','ATR_roi2_R.nii.gz' ),
       'L_CGC':('CGC_roi1_L.nii.gz','CGC_roi2_L.nii.gz' ),
       'R_CGC':('CGC_roi1_R.nii.gz','CGC_roi2_R.nii.gz' ),
       'L_CST':('CST_roi1_L.nii.gz','CST_roi2_L.nii.gz' ),
       'R_CST':('CST_roi1_R.nii.gz','CST_roi2_R.nii.gz' ),
       'FA'   :('FA_L.nii.gz'      ,'FA_R.nii.gz'       ),
       'FP'   :('FP_L.nii.gz'      ,'FP_R.nii.gz'       ),
       'L_HCC':('HCC_roi1_L.nii.gz','HCC_roi2_L.nii.gz' ),
       'R_HCC':('HCC_roi1_R.nii.gz','HCC_roi2_R.nii.gz' ),
       'L_IFO':('IFO_roi1_L.nii.gz','IFO_roi2_L.nii.gz' ),
       'R_IFO':('IFO_roi1_R.nii.gz','IFO_roi2_R.nii.gz' ),
       'L_ILF':('ILF_roi1_L.nii.gz','ILF_roi2_L.nii.gz' ),
       'R_ILF':('ILF_roi1_R.nii.gz','ILF_roi2_R.nii.gz' ),
       'L_SLF':('SLF_roi1_L.nii.gz','SLF_roi2_L.nii.gz' ),
       'R_SLF':('SLF_roi1_R.nii.gz','SLF_roi2_R.nii.gz' ),
       'L_AF' :('SLF_roi1_L.nii.gz','SLFt_roi2_L.nii.gz'),
       'R_AF' :('SLF_roi1_R.nii.gz','SLFt_roi2_R.nii.gz'),
       'L_UNC':('UNC_roi1_L.nii.gz','UNC_roi2_L.nii.gz' ),
       'R_UNC':('UNC_roi1_R.nii.gz','UNC_roi2_R.nii.gz' )}


def makeSureDir(folder):
    if not os.path.exists(folder): os.mkdir(folder)

# Read all the subjects

subs = glob.glob(join(SUBJDIR,'ROIs/s*'))
for i, s in enumerate(subs):
    sname = s.split('/')[-1]
    subs[i] = sname
    
structural = join(SUBJDIR, "Structural")
    
tractsoutdir = join(SUBJDIR, "tracts")

ROIdir = join(SUBJDIR, "ROIs")

makeSureDir(tractsoutdir)

for s in subs:
    makeSureDir(join(tractsoutdir, s))
    makeSureDir(join(tractsoutdir, s, 'tracts-'+fmm))
    makeSureDir(join(tractsoutdir, s, 'tracts-'+fmm,'det'))
    makeSureDir(join(tractsoutdir, s, 'tracts-'+fmm,'prob'))

tcks = join(SUBJDIR, "TCK-Streamlines-"+fmm+"mm")



# Do for all subjects
for s in subs:

    
    # Make sure folders exist
    #if not os.path.exists(join(s,'tracts-'+fmm)): os.mkdir(join(s,'tracts-'+fmm))
    #if not os.path.exists(join(s,'tracts-'+fmm,'det')): os.mkdir(join(s,'tracts-'+fmm,'det'))
    #if not os.path.exists(join(s,'tracts-'+fmm,'prob')): os.mkdir(join(s,'tracts-'+fmm,'prob'))
    
    # See if aparcaseg is binarized
    aparc    = join(structural, s,'T1w','aparc+aseg.nii.gz ')
    binaparc = join(structural, s,'T1w', 'segmentation.nii.gz')
        
    cmd      = str('mrthreshold  -quiet -abs 1000 '+aparc+binaparc);
    
    if not os.path.exists(binaparc): spres = sp.call(cmd, shell=True)
    
    # And now dilate it so that we are sure that the fibers arrive
    binaparcdil = join(structural, s, 'T1w', 'segmentation_dilated.nii.gz')
    cmd      = str('maskfilter -npass 3 '+binaparc+' dilate '+binaparcdil)
    if not os.path.exists(binaparcdil): spres = sp.call(cmd, shell=True)
    
    # Select what we want to use
    if useCortex:
        Cortex = ' -include '+binaparc+' ' # Can be ' '
    else:
        Cortex = ' '
        
    # Get path to the det and prob tractograms
    tractogram = {'det':join(tcks, s,'tracking-deterministic-'+fmm+'mm.tck'),
                  'prob':join(tcks, s,'tracking-probabilistic-'+fmm+'mm.tck')}
                  
    # Obtain tracts dict and loop them
    tracts = R2T.keys()
    for tr in tracts:        
        # Read required ROIs
        ROI1   = join(ROIdir, s,'ROIs', R2T[tr][0])
        ROI2   = join(ROIdir, s,'ROIs', R2T[tr][1])
                
        # Run the mrtrix code
        # Set options
        exclude   = ' ' # can be ' ' or '-exclude roi.nii.gz '
        maxlength = ' ' # can be ' ' or '-maxlength xx '
        minlength = ' ' # can be ' ' or '-minlength xx ' 
        force     = ' -force ' # can be ' ' or '-force '
        # det
        for dodp in doDetProb:
            tracks_in  = tractogram[dodp]
            if useCortex:
                tracks_out = join(tractsoutdir, s,'tracts-'+fmm,dodp,tr+'_cortex.tck')
            else:
                tracks_out = join(tractsoutdir, s,'tracts-'+fmm,dodp,tr+'.tck')
                
            cmd     = str('tckedit -include ' +ROI1+ ' -include ' +ROI2+ ' '+
                             Cortex + exclude + maxlength + minlength + force +
                             tracks_in +' '+tracks_out)
            spres   = sp.call(cmd, shell=True)
            
            
            
            
