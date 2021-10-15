###MIP installer
###Bioinformatics Group, Qingdao University
###Updated at Oct. 14, 2021
###Updated by Xiaoquan Su, Gongchao Jing
#!/bin/bash

###Check Parallel-Meta environment variable###
if [ $ParallelMETA ]
    then
    echo -e "\n**Parallel-Meta is already installed**"
else
    echo -e "\n**Please install latest version of Parallel-Meta**"
    exit
fi

echo -e "\n**MIP Installation**"

###Build source code for src package###
echo -e "\n**MIP src package build**"
make
echo -e "\n**Build Complete**"

###Plugin installation###
echo -e "\n**Plugin Installation**"
cp bin/PM-parse-mip $ParallelMETA/bin
cp -rf mip_16s $ParallelMETA/databases

##Check database configuration##
Check_db_config=`grep -c "mip_16s" $ParallelMETA/databases/db.config`
if [ "$Check_db_config" = "0" ]
    then
    cp $ParallelMETA/databases/db.config $ParallelMETA/databases/.db.config.bk
    cat $ParallelMETA/databases/.db.config.bk db.config > $ParallelMETA/databases/db.config
fi

echo -e "\n**MIP Installation Complete**"
echo -e "\n**An example dataset with demo script is available in \"example\"**"

