## This code was created to predict the fictional authors of Gone Girl, by Gillian Flynn
# The challenge was to differentiate between Nick and Amy based on common stop words
# The data ('authorship.txt') was created by passing the text of GG through stylometric scripts, 
# which were created by:

#   Eder, M., Rybicki, J. (2011). Stylometry with R. In "Digital Humanities 
#   2011: Conference Abstracts." Stanford University, Stanford, CA, pp. 308-11.
library(randomForest)

setwd("~/Documents/Git/WordprintAuthorPrediction/GoneGirl ")

authorship <- read.delim("./authorship.txt")

# Set seed and divide the data into training and test
set.seed(123)
authorship$randu <- runif(63, 0,1)
authorship.train <- authorship[authorship$randu < .4,]
authorship.test <- authorship[authorship$randu >= .4,]

# The model: I tried to avoid all words that would be too predictive prima facie
authorship.model.rf = randomForest(author ~ the + a + to + and + of + it + was + in. + 
                                     that + on + for. + with + is + but + like + be + at + 
                                     so + this + have + what + not + as,
                                   data=authorship.train, ntree=5000, mtry=15, importance=TRUE)

# Make a prediction variable and tables of the outcome
authorship.test$pred.author.rf = predict(authorship.model.rf, authorship.test, type="response")
table(authorship.test$author, authorship.test$pred.author.rf)
prop.table(table(authorship.test$author, authorship.test$pred.author.rf),1)