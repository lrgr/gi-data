echo outputs > processed_outputs.txt
ls datasets/*/output/processed/* >> processed_outputs.txt
chacon create -i processed_outputs.txt -o processed_file_md5_checksums.yml
