echo "######## QIIME pipeline #########"
echo "USAGE:    " $0 " config_file [snakemake options]"
echo "examples: " $0 " config_file"
echo "            Run qiime pipeline with defaults"
echo "          " $0 " config_file -n -p"
echo "            Dont run pipeline; just print out shell command"
echo "          " $0 " configfile deliver"
echo "            deliver the result"
echo "          " $0 " -h"
echo "            Print pipeline options"
echo "          " $0 " "
echo "            "

progPath=`dirname $0`
snakemake_com="/mnt/software/anaconda/envs/py3/bin/snakemake"
snakemake_file=qiime.sm
if hash $snakemake_com 2>/dev/null; then
    if [ -f $1 ]; then
        cfg=$1
        shift
        if [[ " $* " == *" -j "* ]]; then
            set -x
            $snakemake_com -s ${progPath}/${snakemake_file} --configfile ${cfg} "$@"
        else
            let p=`nproc`/2  # use half of the available cores
            set -x
            $snakemake_com -s ${progPath}/${snakemake_file} --configfile ${cfg} "$@" -j ${p}
        fi
    else
        echo "Config file $1 doesn't exist!"
    fi
else
    echo "Cannot find snakemake to run the pipeline"
fi
