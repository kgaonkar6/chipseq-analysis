#!/bash/bin
declare -a array=("Aldoa" "Angpt1" "Bcl2" "Cdkn1a" "Edn1" "Egf" "Eno1" "Epo" "Flt1" "Gapdh" "Hk1" "Hmox1" "Itgb2" "Ldha" "Nos2" "Nos3" "Nppa" "Pdk1" "Pfkfb3" "Pfkl" "Pgk1" "Serpine1" "Slc2a1" "Tek" "Tfrc" "Tie1" "Timp1" "Vegfa")
for i in "${array[@]}";do sed -i '/'"${i}"'/d' $1;done


