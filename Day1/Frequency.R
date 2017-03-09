### --------------------------------------------------------------
### Cat paw example, exact binomial test, pp. 30–31
### http://www.biostathandbook.com/exactgof.html
### --------------------------------------------------------------
### In this example:
###   2 is the number of successes: number of times the cat uses his right paw
###   10 is the number of trials (8 times he uses the left paw)
### can you conclude that he is right-pawed, or could this result have occurred
### due to chance under the null hypothesis that he bats equally with each paw?
###   0.5 is the hypothesized probability of success (unbiased)

dbinom(2, 10, 0.5)            # Probability of single event only!
#   Not binomial test!


binom.test(2, 10, 0.5, 
           alternative="less",       # One-sided test
           conf.level=0.95)

### --------------------------------------------------------------
### Probability density plot, binomial distribution, p. 31
### --------------------------------------------------------------
# In this example:
#   You can change the values for trials and prob
#   You can change the values for xlab and ylab

trials = 10
prob = 0.5

x = seq(0, trials)                   # x is a sequence, 1 to trials
y = dbinom(x, size=trials, p=prob)   # y is the vector of heights

barplot (height=y, 
         names.arg=x,
         xlab="Number of uses of right paw",
         ylab="Probability under null hypothesis")



binom.test(2, 10, 0.5, 
           alternative="two.sided",  # Two-sided test
           conf.level=0.95)



### --------------------------------------------------------------
### Probability density plot, binomial distribution, p. 31
### --------------------------------------------------------------
# In this example:
#   You can change the values for trials and prob
#   You can change the values for xlab and ylab

trials = 10
prob = 0.5

x = seq(0, trials)                   # x is a sequence, 1 to trials
y = dbinom(x, size=trials, p=prob)   # y is the vector of heights

barplot (height=y, 
         names.arg=x,
         xlab="Number of uses of right paw",
         ylab="Probability under null hypothesis")

### --------------------------------------------------------------
### Drosophila example, Chi-square goodness-of-fit, p. 46
### --------------------------------------------------------------
#####
####3:1 ratio of smooth wings to wrinkled wings in offspring from a bunch of Drosophila crosses.
#You observe 770 flies with smooth wings and 230 flies with wrinkled wings; 
#the expected values are 750 smooth-winged and 250 wrinkled-winged flies
#http://www.biostathandbook.com/chigof.html
observed = c(770, 230)        # observed frequencies
expected = c(0.75, 0.25)      # expected proportions

chisq.test(
  x = observed,
  p = expected, 
)



### --------------------------------------------------------------
### Vaccination example, Chi-square independence, pp. 59–60
###      Example directly reading a table as a matrix
### --------------------------------------------------------------

Input =(
  
  "Injection.area  No.severe  Severe
  Thigh           4788       30
  Arm             8916       76
  ")

Matriz = as.matrix(read.table(textConnection(Input),
                              header=TRUE, 
                              row.names=1))

Matriz  




chisq.test(Matriz,
           correct=TRUE)      # Continuity correction for 2 x 2
#      table


chisq.test(Matriz,
           correct=FALSE)      # No continuity correction for 2 x 2
#      table

### --------------------------------------------------------------
### Post-hoc example, Fisher’s exact test, p. 79
### --------------------------------------------------------------
###http://www.biostathandbook.com/fishers.html

Input =(
  
  "Frequency  Damaged  Undamaged
 Daily       1        24
 Weekly      5        20
 Monthly    14        11
 Quarterly  11        14
")

Matriz = as.matrix(read.table(textConnection(Input),
                              header=TRUE, 
                              row.names=1))

Matriz

fisher.test(Matriz,
            alternative="two.sided")
library(RVAideMemoire)
fisher.multcomp(Matriz, 
                p.method = "none")

# Can adjust p-values; 
# See ?p.adjust for options      


