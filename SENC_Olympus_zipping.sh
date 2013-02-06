for i in $(seq 2925 3104)
       do
          sbatch SENC_Olympus_Batch_zipping_unix.txt $i

done

