
basedir=NinN_HCC_GISTIC
segfile=NinN_HCC.seg
refgenefile=/program/install/GISTIC_2_0_23/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat
cnvfile=cnvfile_hg19.refined.txt  # from dgv and nsv database

/program/install/GISTIC_2_0_23/gp_gistic2_from_seg -b $basedir -seg $segfile -refgene $refgenefile -cnv $cnvfile -genegistic 1 -smallmem 1 -broad 1 -armpeel 1 -savegene 1 -ta 0.2 -td 0.2 -qvt 0.25 -cap 1.5 -conf 0.9 -brlen 0.8
