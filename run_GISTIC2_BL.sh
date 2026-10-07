#
##download GISTIC file to the working directory
#wget ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTIC_2_0_23.tar.gz
#tar -xzf GISTIC_2_0_23.tar.gz


#mkdir -p ~/mylibs && ln -sf $(find /usr /lib /lib64 -name "libncurses.so.6" 2>/dev/null | head -1) ~/mylibs/libncurses.so.5 && export LD_LIBRARY_PATH=/main_dir/test_out/cnv_aura/GISTIC_2_0_23/MATLAB_Compiler_Runtime/v83/runtime/glnxa64:/main_dir/test_out/cnv_aura/GISTIC_2_0_23/MATLAB_Compiler_Runtime/v83/bin/glnxa64:/main_dir/test_out/cnv_aura/GISTIC_2_0_23/MATLAB_Compiler_Runtime/v83/sys/os/glnxa64:~/mylibs && GISTIC_DIR=/main_dir/test_out/cnv_aura/GISTIC_2_0_23 && $GISTIC_DIR/gistic2 -b /main_dir/test_out/cnv_aura/gistic_output -seg /main_dir/test_out/cnv_aura/all_samples_gistic.seg -refgene $GISTIC_DIR/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat -genegistic 1 -smallmem 1 -broad 1 -brlen 0.98 -conf 0.95 -armpeel 1 -savegene 1 -gcm extreme -ta 0.1 -td 0.1 -js 15 -rx 0 -v 30
#!/bin/bash

# 
# GISTIC 2 - Full setup and run script
# 

#!/bin/bash
BASE=<your_base_dir>
GISTIC_DIR=$BASE/GISTIC_2_0_23
MCR_DIR=$GISTIC_DIR/MATLAB_Compiler_Runtime/v83
REFGENE=$GISTIC_DIR/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat

IN_DIR=/$BASE/gistic_input_BL
OUT_DIR=/$BASE/gistic_output_BL
mkdir -p $OUT_DIR

export LD_LIBRARY_PATH=$MCR_DIR/runtime/glnxa64:$MCR_DIR/bin/glnxa64:$MCR_DIR/sys/os/glnxa64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/ncurses5_compat:$LD_LIBRARY_PATH
export XAPPLRESDIR=$MCR_DIR/X11/app-defaults

cd $GISTIC_DIR

./gistic2 \
  -b $OUT_DIR \
  -seg $IN_DIR/gistic_segs_BL.txt \
  -mk $IN_DIR/gistic_markers.txt \
  -refgene $REFGENE \
  -genegistic 1 \
  -smallmem 1 \
  -broad 1 \
  -brlen 0.5 \
  -conf 0.90 \
  -qvt 0.10 \
  -armpeel 1 \
  -savegene 1 \
  -gcm extreme \
  -ta 0.1 -td 0.1 \
  -rx 0



