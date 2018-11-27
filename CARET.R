# load random forest,mlbennch dataset, caret
library(randomForest)
library(mlbench)
library(caret)

# load Dataset
data(Sonar)
dataset <- Sonar
x <- dataset[,1:60]
y <- dataset[,61]

# using carets trainControl and train to tune random forest model

# Create model with default paramters
control <- trainControl(method="repeatedcv", number=5, repeats=5)
seed <- 7
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(ncol(x)) # 7
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(Class~., data=dataset, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default)

# Random Search (Tuning Using caret)
control <- trainControl(method="repeatedcv", number=5, repeats=5, search="random")
set.seed(seed)
mtry <- sqrt(ncol(x))
rf_random <- train(Class~., data=dataset, method="rf", metric=metric, tuneLength=15, trControl=control)
print(rf_random)
plot(rf_random)
