if [ -z "$1" ]
      then
          echo "No release argument supplied"
          echo "USAGE: sh <script> <release>"
          exit
fi


target=lrgr:gi-data/$1
echo "Deleting $target for clean upload"

rmobj -r $target

for d in \
collins_et_al \
roguev_et_al
do
    echo "Copying $d to ObjStore"
    cpobj -r all_processed_datasets/$d lrgr:gi-data/$1
done


echo "Uploading metadata (git commit etc.) to ObjStore"
cpobj all_processed_datasets/meta.txt lrgr:gi-data/$1/meta.txt
cpobj all_processed_datasets/checksums.yml lrgr:gi-data/$1/checksums.txt
cpobj all_processed_datasets/verify_result.yml lrgr:gi-data/$1/verify_result.txt
