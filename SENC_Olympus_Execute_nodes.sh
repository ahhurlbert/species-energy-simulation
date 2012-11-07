for i in $(seq 1 2592)
       do
          sbatch SENC_Olympus_Batch_nodes_unix.txt $i

done

