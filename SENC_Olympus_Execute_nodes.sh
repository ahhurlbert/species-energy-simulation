for i in $(seq 2927 3104)
       do
          sbatch SENC_Olympus_Batch_nodes_unix.txt $i

done

