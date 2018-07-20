mastercommand <- sapply(1:100, function(i){
  print(paste0("sbatch submitscript", i, ".txt;"))
})
cat(mastercommand)
cat(mastercommand, file = "masterscript.txt")

