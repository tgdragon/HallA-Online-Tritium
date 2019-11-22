source /apps/root/6.18.00/setroot_CUE
setenv ANALYZER /work/halla/triton/ole/analyzer
setenv LD_LIBRARY_PATH ${ANALYZER}/lib64:${LD_LIBRARY_PATH}
setenv PATH ${ANALYZER}/bin::${PATH}
setenv DB_DIR `pwd`/DB
