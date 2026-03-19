FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04 
#AS builder

ENV CUDA_HOME=/usr/local/cuda-11.8
ENV PATH=/usr/local/cuda-11.8/bin:$PATH
ENV LD_LIBRARY_PATH=/usr/local/cuda-11.8/lib64:$LD_LIBRARY_PATH

RUN  apt update && apt install -y \
        git wget curl \
        build-essential \
        cmake \
        libncurses5-dev \
        zlib1g-dev liblzma-dev libbz2-dev libcurl4-openssl-dev libspdlog-dev libboost-all-dev -y \
        && rm -rf /var/lib/apt/list

# -----------------------------------------------------------------
# Install micromamba
# -----------------------------------------------------------------

RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba; 
RUN mkdir -p /opt/micromamba
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
RUN /bin/micromamba shell init -s bash -r /opt/micromamba

# -----------------------------------------------------------------
# Install general python env
# -----------------------------------------------------------------
RUN echo 'name: benchmark_env\n\
channels:\n\
- conda-forge\n\
- anaconda\n\
- bioconda\n\
- defaults\n\
- https://repo.anaconda.com/pkgs/main\n\
- https://repo.anaconda.com/pkgs/r\n\
dependencies:\n\
- biopython=1.86\n\
- python\n\
- pip\n\
- polars=1.36.1\n\
- seaborn=0.13.2\n\
- snakemake-minimal=9.14.5\n\
- samtools==1.21\n\
- bedtools==2.30.0\n\
- nanoplot==1.46.2\n\
- nanoq==0.10.0\n\
- nanostat==1.6.0\n\
- nanoplot=1.46.2\n\
- minimap2==2.28.0\n\
- time==1.9.0\n\
- pip:\n\
    - blue-crab==0.4.0' > benchmark_env.yml; \
    micromamba env create -f benchmark_env.yml; \
    echo 'micromamba activate benchmark_env' >> /root/.bashrc; \
    rm  benchmark_env.yml

# -----------------------------------------------------------------
# HSTlib
# -----------------------------------------------------------------

ADD ./docker_build_assets/htslib-1.21.tar.bz2 /tooling
RUN cd /tooling/htslib-1.21; \
    ./configure --prefix=/usr/local; \
    make -j && make install; 

# -----------------------------------------------------------------
# DeepBAM
# -----------------------------------------------------------------

## DeepBAM conda setup
ADD docker_build_assets/DeepBAM.tar.gz /tooling
RUN eval $(micromamba shell hook --shell bash) && \
    micromamba create -n DeepBAM python=3.11 -y && \
    micromamba run -n DeepBAM pip install numpy torch==2.0.1

## DeepBAM build
SHELL ["/bin/micromamba", "run", "-n", "DeepBAM", "/bin/bash", "-c"]
RUN mkdir -p /tooling/DeepBAM/cpp/build && \
    cd /tooling/DeepBAM/cpp/build && \
    export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH && \
    export CMAKE_PREFIX_PATH=`python -c 'import torch;print(torch.utils.cmake_prefix_path)'` && \
    cmake -DCMAKE_PREFIX_PATH=`python -c 'import torch;print(torch.utils.cmake_prefix_path)'` .. && \
    make -j

## setup DeepBAM executable
SHELL ["/bin/bash", "-c"]
RUN echo -e '#!/bin/bash\n\
eval "$(micromamba shell hook --shell bash)"\n\
micromamba activate DeepBAM\n\
export CUDA_HOME="/usr/local/cuda-11.8"\n\
export LD_LIBRARY_PATH="/usr/local/cuda-11.8/lib64:$LD_LIBRARY_PATH"\n\
export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"\n\
/tooling/DeepBAM/cpp/build/DeepBAM $@;' > /usr/local/bin/DeepBAM && \
    chmod +x /usr/local/bin/DeepBAM;   

# -----------------------------------------------------------------
# DeepPlant
# -----------------------------------------------------------------

## DeepPlant conda setup

ADD docker_build_assets/DeepPlant.tar.gz /tooling

RUN eval $(micromamba shell hook --shell bash) && \
    micromamba create -n DeepPlant python=3.12 -y && \
    micromamba run -n DeepPlant pip install numpy==1.26 torch==2.2.1

SHELL ["/bin/micromamba", "run", "-n", "DeepBAM", "/bin/bash", "-c"]
RUN mkdir /tooling/DeepPlant/build && cd $_ && \
    cmake -DCMAKE_PREFIX_PATH=`python -c 'import torch;print(torch.utils.cmake_prefix_path)'` ..  && \
    export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH && \
    make -j

## setup DeepPlant executable
SHELL ["/bin/bash", "-c"]
RUN echo -e '#!/bin/bash\n\
eval "$(micromamba shell hook --shell bash)"\n\
micromamba activate DeepPlant\n\
export CUDA_HOME="/usr/local/cuda-11.8"\n\
export LD_LIBRARY_PATH="/usr/local/cuda-11.8/lib64:$LD_LIBRARY_PATH"\n\
export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"\n\
/tooling/DeepPlant/build/DeepPlant $@;' > /usr/local/bin/DeepPlant && \
    chmod +x /usr/local/bin/DeepPlant;  

# -----------------------------------------------------------------
# Rockfish
# -----------------------------------------------------------------
ADD docker_build_assets/rockfish.tar.gz /tooling
RUN eval $(micromamba shell hook --shell bash) && \
    cd /tooling/rockfish && \
    micromamba create -n rockfish python=3.9  && \
    micromamba run -n rockfish pip install --extra-index-url https://download.pytorch.org/whl/cu118 . && \
    echo -e '#!/bin/bash\n\
eval "$(micromamba shell hook --shell bash)"\n\
micromamba activate rockfish\n\
rockfish $@' > /usr/local/bin/rockfish && \
    chmod +x /usr/local/bin/rockfish 

# -----------------------------------------------------------------
# f5c
# -----------------------------------------------------------------   

COPY docker_build_assets/f5c-v1.5-binaries.tar.gz /tooling
RUN mkdir /tooling/f5c-v1.5 && cd /tooling/f5c-v1.5/ && \
    tar xvf ../f5c-v1.5-binaries.tar.gz -C . --strip-component=1  && \
    rm ../f5c-v1.5-binaries.tar.gz &&\
    echo -e '#!/bin/bash\n\
    /tooling/f5c-v1.5/f5c_x86_64_linux_cuda $@' > /usr/local/bin/f5c_x86_64_linux_cuda && \
    chmod +x /usr/local/bin/f5c_x86_64_linux_cuda && \
    echo -e '#!/bin/bash\n\
    /tooling/f5c-v1.5/f5c_x86_64_linux $@' > /usr/local/bin/f5c_x86_64_linux && \
    chmod +x /usr/local/bin/f5c_x86_64_linux && \
    echo -e '#!/bin/bash\n\
    /tooling/f5c-v1.5/f5c_x86_64_linux $@' > /usr/local/bin/f5c && \
    chmod +x /usr/local/bin/f5c;

# -----------------------------------------------------------------
# Dorado
# ----------------------------------------------------------------- 

## dorado install
RUN mkdir /tooling/dorado -p
ADD ./docker_build_assets/dorado-0.9.1-linux-x64.tar.gz /tooling/dorado
ADD ./docker_build_assets/dorado-1.1.1-linux-x64.tar.gz /tooling/dorado

## download dorad-0.9.1
RUN cd /tooling/dorado; \
    echo -e '#!/bin/bash\n\
    /tooling/dorado/dorado-0.9.1-linux-x64/bin/dorado $@' > /usr/local/bin/dorado-0.9.1 && \
    chmod +x /usr/local/bin/dorado-0.9.1

## download dorad-1.1.1
RUN echo -e '#!/bin/bash\n\
    /tooling/dorado/dorado-1.1.1-linux-x64/bin/dorado $@' > /usr/local/bin/dorado-1.1.1 && \
    chmod +x /usr/local/bin/dorado-1.1.1

# -----------------------------------------------------------------
# Modkit
# -----------------------------------------------------------------
ADD ./docker_build_assets/modkit_v0.5.1rc1_u16_x86_64.tar.gz /tooling/modkit
RUN cd /tooling/modkit/ && \
    mv dist_modkit_v0.5.1_8fa79e3/* ./ && \
    rm -r dist_modkit_v0.5.1_8fa79e3 && \
    echo -e '#!/bin/bash\n\
    /tooling/modkit/modkit $@' > /usr/local/bin/modkit && \
    chmod +x /usr/local/bin/modkit

# # -----------------------------------------------------------------
# # Rerio
# # -----------------------------------------------------------------

ADD ./docker_build_assets/rerio.tar.gz /tooling
RUN cd /tooling/rerio && \
    rm -r basecall_models clair3_models remora_models taiyaki_models ONT_logo.png README.rst; \
    echo -e '#!/bin/bash\n\
    /tooling/rerio/download_model.py --dorado; \n\
    cp -r /tooling/rerio/dorado_models/*/ $@' > /usr/local/bin/download_model.py; \
    chmod +x /usr/local/bin/download_model.py

# -----------------------------------------------------------------
# DeepMod2
# -----------------------------------------------------------------    
ADD ./docker_build_assets/deepmod2.tar.gz /tooling
RUN cd /tooling/deepmod2  && \
    micromamba create -y -n deepmod2 \
        -c pytorch -c nvidia -c conda-forge \
        pytorch torchvision pytorch-cuda=11.8 && \
    micromamba env update -n deepmod2 -f environment.yml -y

RUN echo '#!/bin/bash\n\
eval "$(micromamba shell hook --shell bash)"\n\
micromamba activate deepmod2\n\
python /tooling/deepmod2/deepmod2 $@' > /usr/local/bin/deepmod2 && \
    chmod +x /usr/local/bin/deepmod2;    

# ----------------------------------------------------------------

## cleanup
RUN micromamba clean --all --yes
RUN rm -r /tooling/rockfish
RUN rm -r /tooling/DeepPlant/model/
RUN rm -r /tooling/DeepBAM/traced_script_module/

# # Stage 2: Runtime stage - minimal!
# FROM nvidia/cuda:11.8.0-base-ubuntu22.04

# RUN mkdir /tooling
# # # Copy only the Python packages from builder
# COPY --from=builder /tooling /
# COPY --from=builder /opt/micromamba /opt
# COPY --from=builder /usr/local/bin/* /usr/local/bin/*
# # COPY --from=builder /usr/local/bin /usr/local/bin


ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV MAMBA_EXE=/bin/micromamba
ENV CONDA_PREFIX=/opt/micromamba/envs/benchmark_env
ENV CONDA_DEFAULT_ENV=benchmark_env
ENV PATH=/opt/micromamba/envs/benchmark_env/bin:/opt/micromamba/bin:$PATH
ENV CONDA_SHLVL=1

# RUN echo -e '#!/bin/bash\n\
# figlet "Sowpati LAB";\n\
# echo "-------------------------------------------------------------";\n\
# echo "              ont-methylation-benchmark-env                  ";\n\
# echo "-------------------------------------------------------------";\n\
# echo "" \n\
# "$@"' > /usr/local/bin/entrypoint; \
#     chmod +x /usr/local/bin/entrypoint
# eval "$(micromamba shell hook --shell bash)"; \n\
# micromamba activate benchmark_env; \n\

# ENTRYPOINT ["entrypoint"]
CMD ["/bin/bash"]