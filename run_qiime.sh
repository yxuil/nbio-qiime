echo "######## QIIME pipeline #########"
echo "USAGE:    " $0 " [snakemake options]"
echo "examples: " $0 " "
echo "            Run qiime pipeline with defaults"
echo "          " $0 " -n -p"
echo "            Dont run pipeline; just print out shell command"
echo "          " $0 " --configfile [configfile]"
echo "            Specify run config file to use"
echo "          " $0 " -h"
echo "            Print pipeline options"
#echo "          " $0 " "
#echo "            "

progPath=`dirname $0`
snakemake_com="/mnt/software/anaconda/envs/py3/bin/snakemake"
snakemake_file=qiime.sm
if hash $snakemake_com 2>/dev/null; then
    set -x
    let p=`nproc`/2
    $snakemake_com -s ${progPath}/${snakemake_file} "$@" -j ${p}
else
    echo "Cannot find snakemake to run the pipeline"
fi
