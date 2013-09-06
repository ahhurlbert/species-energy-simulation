for i in $(seq 4065 4067)
       do
          sbatch SENC_Olympus_Batch_analysis_unix.txt $i
done

