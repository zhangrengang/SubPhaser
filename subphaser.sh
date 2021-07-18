#/bin/bash
config=sg.config
genome=$1
shift
opts=$@
python3 /share/home/nature/users/zhangrenang/subphaser/subphaser.py -i $genome -c $config $opts
