
d <- read.table("MYPROFLEFILE.csv", header=T, sep=";");
plot( d[,1], d[,3], xlab="sequence index", ylab="minimal energy", type="l", col="blue", lwd=2)
abline(h=0, col="red", lty=2, lwd=2)


# read data, skip first column, and replace NA and E>0 values with 0
d <- read.table("MYPAIRMINE.csv",header=T,sep=";");
d <- d[,2:ncol(d)];
d[is.na(d)] = 0;
d[d>0] = 0;
# plot
image( 1:nrow(d), 1:ncol(d), as.matrix(d), col = heat.colors(100), xlab="index in sequence 1", ylab="index in sequence 2");
box();



library(tidyverse)

read_delim("pairMinE-t1q1.csv", delim=";") %>%
  replace(is.na(.), 0) %>%
  rename( target = minE ) %>%
  pivot_longer( cols=-target, names_to = "query", values_to = "E" ) %>%
  separate( query, into = c("qSeq","qIdx"), sep = "_") %>%
  separate( target, into = c("tSeq","tIdx"), sep = "_") %>%
  mutate(across(c(tIdx,qIdx,E),as.numeric)) %>%
  mutate(E = ifelse( E>0, 0, E)) %>%
  ggplot(aes(x=tIdx,y=qIdx, fill=E)) +
  geom_tile() + 
  scale_fill_gradient(low="blue",high="red")