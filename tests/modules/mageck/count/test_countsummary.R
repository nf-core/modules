Sweave("test_countsummary.Rnw");
library(tools);

texi2dvi("test_countsummary.tex",pdf=TRUE);

