# Wrapper script to submit QIIME pipeline to SGE cluster
# Xiaoyu Liu
# 2014.1-2015.1

USAGE="Usage: run_qiime_sge.sh run_config_file [-c \"qsub_command\"] [snakemake options]
  run_config_file: Run configuration in JSON format. Refer to the example.
                   Required paramter
  qsub_command:    qsub command, like: \"qsub -V -cwd -j y -b y -o {log} -q all.q\"
                   if it is specified, it has to be right after config file
  options:         Other valid snakemake options"
# if no qsub commands specified by -c option, use default:
qsub_default="qsub -V -cwd -j y -b y -o {log} -q all.q"

if hash qsub 2>/dev/null; then
    # check snakemake installation. It should be updated whenever snakemake is
    # changed
    progPath=`dirname $0`
    snakemake_com=/mnt/software/anaconda/envs/py3/bin/snakemake
    if hash $snakemake_com 2>/dev/null; then
        if [ $# -eq 0 ]; then  # no config file specified
            echo "$USAGE"
        else
            if [ "$1" == "-h" ]; then  # print help message
                echo  "$USAGE"
                echo
                echo "Options for snakemake:\n"
                $snakemake_com -h
            elif [ -f $1 ]; then # config file exist, but no grid submission command
                config=$1
                shift
                set -x
                $snakemake_com --cluster "${qsub_default}" --jn s.\{rulename\}.\{jobid\}.sh -s ${progPath}/qiime.sm \
                --configfile ${config} "$@" -j 10 -w 30
            else
                echo "config file $1 does not exist!"
                echo "$USAGE"
            fi
        fi
    else
        echo "Cannot find snakemake at " $snakemake_com "!"
    fi
else
    echo "Cannot find command qsub. Use run_qiime.sh instead"
fi
