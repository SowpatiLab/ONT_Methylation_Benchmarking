[ ! -d docker_build_assets ] && mkdir docker_build_assets;
cd docker_build_assets

[ ! -f dorado-0.9.1-linux-x64.tar.gz ] && wget  https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.1-linux-x64.tar.gz
[ ! -f dorado-1.1.1-linux-x64.tar.gz ] && wget  https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.1.1-linux-x64.tar.gz
[ ! -f modkit_v0.5.1rc1_u16_x86_64.tar.gz ] && wget https://github.com/nanoporetech/modkit/releases/download/v0.5.1-rc1/modkit_v0.5.1rc1_u16_x86_64.tar.gz
if [ ! -f deepmod2.tar.gz ];
then 
    git clone https://github.com/WGLab/DeepMod2.git deepmod2;
    rm -rf deepmod2/.git
    cat deepmod2/src/detect.py | perl -pe 's/^(from ont_fast5_api.*)/#$1/' > tmp.txt
    cat tmp.txt > deepmod2/src/detect.py
    tar -cvzf deepmod2.tar.gz deepmod2
    rm -rf deepmod2 tmp.txt
fi

if [ ! -f rerio.tar.gz ]; 
then
    git clone https://github.com/nanoporetech/rerio.git
    rm -rf rerio/.git
    tar -cvzf  rerio.tar.gz rerio
    rm -rf rerio
fi

if [ ! -f rockfish.tar.gz ]; 
then
    git clone -b r10.4.1 https://github.com/lbcb-sci/rockfish.git --single-branch
    rm -rf rockfish/.git
    tar -cvzf  rockfish.tar.gz rockfish
    rm -rf rockfish
fi

if [ ! -f DeepBAM.tar.gz ]; 
then
    git clone https://github.com/xiaochuanle/DeepBAM.git
    rm -rf DeepBAM/.git
    tar -cvzf  DeepBAM.tar.gz DeepBAM
    rm -rf DeepBAM
fi

if [ ! -f DeepPlant.tar.gz ]; 
then
    git clone https://github.com/xiaochuanle/DeepPlant.git
    rm -rf DeepPlant/.git
    tar -cvzf  DeepPlant.tar.gz DeepPlant
    rm -rf DeepPlant
fi

[ ! -f f5c-v1.5-binaries.tar.gz ] && wget https://github.com/hasindu2008/f5c/releases/download/v1.5/f5c-v1.5-binaries.tar.gz
[ ! -f miniconda.sh ] && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
[ ! -f htslib-1.21.tar.bz2 ] && wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2

cd ../
docker build --no-cache -t sowpati/ont-methylation-benchmarking:latest .
