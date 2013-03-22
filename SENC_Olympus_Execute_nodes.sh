for i in $(seq 3405 3444)
       do
          sbatch SENC_Olympus_Batch_nodes_unix.txt $i
	sleep 1
done

