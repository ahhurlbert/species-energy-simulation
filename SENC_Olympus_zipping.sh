for i in $(1)
       do
          sbatch SENC_Olympus_Batch_zipping_unix.txt $i

done

