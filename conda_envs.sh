git clone git@github.com:brandonlind/pythonimports.git
ln -s pythonimports $HOME/pythonimports

unset PYTHONPATH
unset PYTHONHOME

conda config --add channels conda-forge
conda config --add channels bioconda

## pop_struct env - python pipeline
conda create -n pop_struct python=3.12.2 --yes

conda run --no-capture-output -n pop_struct pip install numpy pandas matplotlib session_info matplotlib_venn ipyparallel cartopy tqdm paramiko pdf2image seaborn geopandas notebook jupyterlab_widgets ipywidgets

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

conda run --no-capture-output -n lea_bigsnpr conda install -y bioconductor-lea=3.6.0

conda run --no-capture-output -n lea_bigsnpr Rscript -e 'install.packages("bigsnpr", repos=c("https://privefl.r-universe.dev","https://cloud.r-project.org"))'

## bcftools 1.23
conda create -n bcftools bioconda::bcftools --yes

## lostruct
conda create -n lostruct -c conda-forge r-base r-devtools gcc_linux-64 gxx_linux-64 make zlib bzip2 xz liblzma-devel libcurl --yes

conda run --no-capture-output -n lostruct R -e "
  install.packages('data.table', repos='https://cloud.r-project.org');
  devtools::install_github('petrelharp/local_pca/lostruct', upgrade='always');
"

## hierfstat
conda create -n hierfstat -c conda-forge r-hierfstat --yes

