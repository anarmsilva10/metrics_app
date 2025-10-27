FROM python:3.12.2-slim

# Set working directory
WORKDIR /metrics_app

# Instal dependecies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda
RUN curl -L https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda.sh&&\
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    /opt/conda/bin/conda clean --all

# Add Conda to the PATH
ENV PATH=/opt/conda/bin:$PATH

# Copy project files into the container
COPY . .

# Accept TOS for Conda channels
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Create and activate the conda environment
RUN conda env create -f environment.yml && conda clean -a

# Expose Port
EXPOSE 8000

# Entrypoint
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "quality_app_env", "shiny", "run", "app.py", "--host=0.0.0.0", "--port=8000"]