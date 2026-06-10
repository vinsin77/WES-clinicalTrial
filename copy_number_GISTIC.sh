#
##download GISTIC file to the working directory
#wget ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTIC_2_0_23.tar.gz
#tar -xzf GISTIC_2_0_23.tar.gz


#mkdir -p ~/mylibs && ln -sf $(find /usr /lib /lib64 -name "libncurses.so.6" 2>/dev/null | head -1) ~/mylibs/libncurses.so.5 && export LD_LIBRARY_PATH=/main_dir/test_out/cnv_aura/GISTIC_2_0_23/MATLAB_Compiler_Runtime/v83/runtime/glnxa64:/main_dir/test_out/cnv_aura/GISTIC_2_0_23/MATLAB_Compiler_Runtime/v83/bin/glnxa64:/main_dir/test_out/cnv_aura/GISTIC_2_0_23/MATLAB_Compiler_Runtime/v83/sys/os/glnxa64:~/mylibs && GISTIC_DIR=/main_dir/test_out/cnv_aura/GISTIC_2_0_23 && $GISTIC_DIR/gistic2 -b /main_dir/test_out/cnv_aura/gistic_output -seg /main_dir/test_out/cnv_aura/all_samples_gistic.seg -refgene $GISTIC_DIR/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat -genegistic 1 -smallmem 1 -broad 1 -brlen 0.98 -conf 0.95 -armpeel 1 -savegene 1 -gcm extreme -ta 0.1 -td 0.1 -js 15 -rx 0 -v 30
#!/bin/bash

# 
# GISTIC 2 - Full setup and run script
# 

BASE=/main_dir/test_out/cnv_aura
GISTIC_DIR=$BASE/GISTIC_2_0_23
MCR_DIR=$GISTIC_DIR/MATLAB_Compiler_Runtime/v83
SEG_FILE=$BASE/all_samples_gistic.seg
REFGENE=$GISTIC_DIR/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat
OUT_DIR=$BASE/gistic_output

# 
# Step 1: Fix libncurses.so.5 
# 
mkdir -p ~/mylibs
NCURSES6=$(find /usr /lib /lib64 -name "libncurses.so.6" 2>/dev/null | head -1)
if [ -z "$NCURSES6" ]; then
    echo "ERROR: libncurses.so.6 not found on system"
    exit 1
fi
ln -sf $NCURSES6 ~/mylibs/libncurses.so.5
echo "libncurses.so.5 symlinked from $NCURSES6"

# 
# Step 2: Set MCR library paths
# 
export LD_LIBRARY_PATH=$MCR_DIR/runtime/glnxa64:$MCR_DIR/bin/glnxa64:$MCR_DIR/sys/os/glnxa64:~/mylibs
export XAPPLRESDIR=$MCR_DIR/X11/app-defaults

# 
# Step 3: Create output directory
# 
mkdir -p $OUT_DIR

# 
# Step 4: Run GISTIC 2
# 
echo "Starting GISTIC 2..."
echo "  Seg file : $SEG_FILE"
echo "  Refgene  : $REFGENE"
echo "  Output   : $OUT_DIR"
echo ""

$GISTIC_DIR/gistic2 \
  -b $OUT_DIR \
  -seg $SEG_FILE \
  -refgene $REFGENE \
  -genegistic 1 \
  -smallmem 1 \
  -broad 1 \
  -brlen 0.98 \
  -conf 0.95 \
  -armpeel 1 \
  -savegene 1 \
  -gcm extreme \
  -ta 0.1 \
  -td 0.1 \
  -js 15 \
  -rx 0 \
  -v 30

# 
# Step 5: Check output
# 
if [ $? -eq 0 ]; then
    echo ""
    echo "GISTIC 2 finished successfully!"
    echo "Output files:"
    ls -lh $OUT_DIR
else
    echo ""
    echo "GISTIC 2 failed — check error messages above"
fi



