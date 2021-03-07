# CDFs-for-Gene-Pathway-Expression
This project's goal was to compare the interferon pathway expression of cells in our disease samples to assess if the disease correlated with over expression in this pathway or under expression. This was based off of a previously published paper's methods which I recreated by the description.

The first few lines are simply pruning the data to the genes and subjects of interest. 

I then performed some quality control on the samples to ensure we were taking reads that were coming through strongly and wouldn't simply be noise. 

The rest of the code is creating the interferon scores, normalizing them, and calculating their cumulative probabilities to be plotted. Those results are included in this folder with the figure description. I also calculated the differences between the two curves and evaluted if the difference was statistically significant. 
