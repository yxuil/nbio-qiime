USAGE="run_qiime_sge.sh [\"qsub_command\"]\nPut qsub command inside quote, like: \"qsub -V -cwd -j y -b y -o {log} -q brchigh.q\""

if hash qsub 2>/dev/null; then
    progPath=`dirname $0`
    snakemake_com=/mnt/software/anaconda/envs/py3/bin/snakemake
    if hash $snakemake_com 2>/dev/null; then
        if [ "$1" == "-h" ]; then
            echo $USAGE
        elif [ $# -gt 2 ]; then
            qsub_opt=$[1]
            shift
        elif [ $# -eq 1 ]; then
            qsub_opt="qsub -V -cwd -j y -b y -o {log} -q brchigh.q"
        fi
        $snakemake_com --cluster "$qsub_opt" --jn s.\{rulename\}.\{jobid\}.sh -s ${progPath}/qiime.sm "$@"

    else
        echo "Cannot find snakemake!"
    fi
else
    echo "Cannot find command qsub. Use run_qiime.sh instead"
fi