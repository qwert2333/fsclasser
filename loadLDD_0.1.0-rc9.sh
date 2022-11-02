unset MARLIN_DLL
#source /cefs/higgs/guofy/cepcsoft/init_ilcsoft.sh
source /cvmfs/cepc.ihep.ac.cn/software/cepcenv/setup.sh
cepcenv use 0.1.1

#Dir_PATH=$(cd "$(dirname ".")";pwd)
export FSClasser_HOME=/cefs/higgs/guofy/CEPC/FSClasser/
export LD_LIBRARY_PATH="$FSClasser_HOME/lib:$LD_LIBRARY_PATH"
export MARLIN_DLL=$MARLIN_DLL:$FastJet_LIB/libfastjet.so:$FSClasser_HOME/lib/libFSClasser.so

