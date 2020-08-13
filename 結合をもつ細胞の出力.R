#Œ‹‡‚ð‚à‚Â×–E‚Ìo—Í
sink("Œ‹‡‚Ì•ûŒü«.txt")
i=1
while (i<=ROI) {
  j=1
  while (j<=ROI) {
    if(FC[i,j]==1){
      cat(i,"¨",j,"\n")
    }
    j=j+1
  }
  i=i+1
}
sink()