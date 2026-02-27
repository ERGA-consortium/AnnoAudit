FROM mambaorg/micromamba
ADD environment.yml /tmp/environment.yml
RUN micromamba create -f /tmp/environment.yml -y && micromamba clean -afy
RUN echo "source activate omark" > ~/.bashrc
ENV PATH /opt/conda/envs/omark/bin:$PATH