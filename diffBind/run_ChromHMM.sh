 JAVA7/java -mx6000M -jar $ChromHMM/ChromHMM.jar BinarizeBam $ChromHMM/CHROMSIZES/${refgenome}.txt $input_dir/mapout $input_dir/cellmarkfile.txt $ChromHMM_binary
JAVA7/java -Xmx10G -jar $ChromHMM/ChromHMM.jar LearnModel -p 1 $ChromHMM_binary $ChromHMM_output $state $refgenom

