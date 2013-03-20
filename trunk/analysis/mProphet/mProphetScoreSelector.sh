#!/bin/bash
[ "$#" -eq 0 ] && echo "Usage: $0 file main_var [vars...]" && exit

file=$1
#remove previous var_ annotations, add main_var to sedstring
sedstring="1s/main_var_//;1s/var_//g;1s/$2/main_var_$2/;"
shift 2
for i in $@
do
    #add vars to sedstring
	sedstring="${sedstring}1s/${i}/var_${i}/;"
done
#execute sed command once only
sed -i"" "$sedstring" $file
