# Use the official Python 3.12 slim image as base
FROM python:3.12-slim


RUN apt-get update && apt-get install -y git-all wget bzip2 && rm -rf /var/lib/apt/lists/*
ARG MINICONDA_VERSION="latest"
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -O /tmp/miniconda.sh

# Install Miniconda silently to /opt/conda
RUN bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Add Conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

RUN conda config --remove channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

COPY environment.yml .  
RUN conda env create -f environment.yml
# RUN conda clean -a -y
# RUN conda init
SHELL ["conda", "run", "-n", "bohra3", "/bin/bash", "-c"]
RUN pip3 install git+https://github.com/MDU-PHL/bohra.git@rethink_structure

RUN bohra check --install-deps --no-databases

# # CMD ["conda", "run", "-n", "bohra3", "bohra", "--help"]