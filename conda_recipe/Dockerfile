FROM continuumio/miniconda
MAINTAINER Remi Lehe <rlehe@lbl.gov>

RUN apt-get update \
    && apt-get install -y \
        gcc \
        libgl1-mesa-glx \
    && rm -rf /var/lib/apt/lists/*

RUN conda install --yes conda conda-build

CMD cd /home/ \
    && conda build --python=2.7 . \
    && conda build --python=3.4 . \
    && conda build --python=3.5 . \
    && conda build --python=3.6 . \
    && conda convert $(conda build --python=2.7 . --output) -p osx-64 \
    && conda convert $(conda build --python=2.7 . --output) -p linux-64 \
    && conda convert $(conda build --python=3.4 . --output) -f -p osx-64 \
    && conda convert $(conda build --python=3.4 . --output) -f -p linux-64 \
    && conda convert $(conda build --python=3.5 . --output) -f -p osx-64 \
    && conda convert $(conda build --python=3.5 . --output) -f -p linux-64 \
    && conda convert $(conda build --python=3.6 . --output) -f -p osx-64 \
    && conda convert $(conda build --python=3.6 . --output) -f -p linux-64
