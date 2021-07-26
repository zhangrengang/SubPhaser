#/bin/bash
config=sg.config
genome=$1
shift
opts=$@
pre=$(basename $(pwd))_
python3 /share/home/nature/users/zhangrenang/subphaser/test_subphaser.py -i $genome -c $config $opts -pre $pre
