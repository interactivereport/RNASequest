#!/usr/bin/env bash
appEnvName="RNAsequest"

condaPath=$(which conda)
if [[ ${#condaPath} -lt 3 ]]; then
    echo "Missing conda"
    echo "Please install conda and add it into PATH"
    exit
else
    echo "conda in $condaPath"
fi

src="$(dirname $0)/src"
condaPath=$(dirname $(dirname $condaPath))
if { conda env list | grep "^$appEnvName"; } >/dev/null 2>/dev/null; then conda env remove -n $appEnvName; fi
#conda env create -f install.yml
# mamba is not in the base conda
conda create -y -n $appEnvName "mamba=0.24.0" -c conda-forge
source $condaPath/etc/profile.d/conda.sh
conda activate $appEnvName
mamba env update -f install.yml

echo "export condaEnv='. $condaPath/etc/profile.d/conda.sh;conda activate $CONDA_PREFIX'" > $src/.env
echo "export PATH=$PATH" >> $src/.env
echo "export OPENBLAS_NUM_THREADS=1" >> $src/.env
echo "export SGE_EXECD_PORT=$SGE_EXECD_PORT" >> $src/.env
echo "export SGE_QMASTER_PORT=$SGE_QMASTER_PORT" >> $src/.env
echo "export SGE_ROOT=$SGE_ROOT" >> $src/.env
echo "export SLURM_CONF=$SLURM_CONF" >> $src/.env
conda deactivate

## additional packages which are not available on anaconda
env -i src="$src" bash -c 'source $src/.env;eval $condaEnv;R -q -e '"'"'suppressWarnings({if(!require(revealjs)) install.packages("revealjs",repos = "https://cran.us.r-project.org")})'"'"

echo "If no errors above, RNAsequest installation/configuration is successful!"
