
### Quick install and start ###
```
git clone https://github.com/zhangrengang/SubPhaser
cd SubPhaser

# install
conda env create -f SubPhaser.yaml
conda activate SubPhaser
python setup.py install

# start
cd example_data
# small genome (Arabidopsis_suecica: 270Mb)
bash test_Arabidopsis.sh
# middle genome	(peanut: 2.6Gb)
bash test_peanut.sh
# large genome (wheat: 14Gb)
bash test_wheat.sh
```
