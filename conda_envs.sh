unset PYTHONPATH
unset PYTHONHOME

conda config --add channels conda-forge
conda config --add channels bioconda

## pop_struct env - python pipeline
conda create -n pop_struct python=3.12.2 --yes
conda activate pop_struct
pip install numpy pandas matplotlib session_info matplotlib_venn ipyparallel cartopy tqdm paramiko pdf2image seaborn geopandas notebook jupyterlab_widgets ipywidgets

## lea_bigsnpr - LEA, bigsnpr, PLINK
conda create -n lea_bigsnpr -c conda-forge \
r-base=4.1.* \
r-essentials \
r-biocmanager \
r-remotes \
r-data.table \
r-rcpp \
libgomp \
plink -y

conda activate lea_bigsnpr
conda install bioconductor-lea=3.6.0

Rscript -e 'install.packages("bigsnpr", repos=c("https://privefl.r-universe.dev","https://cloud.r-project.org"))'

## bcftools
bcftools 1.23
conda create -n bcftools bioconda::bcftools --yes

## lostruct
conda create -n lostruct -c conda-forge r-devtools r-data.table --yes
conda activate lostruct
R -e "devtools::install_github('petrelharp/local_pca/lostruct')"

