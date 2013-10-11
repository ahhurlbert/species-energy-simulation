for i in $(seq 3465 3664)
       do
          sbatch SENC_Olympus_Batch_analysis_unix.txt $i
done

