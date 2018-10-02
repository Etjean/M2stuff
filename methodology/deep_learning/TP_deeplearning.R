library(keras)
install_keras()



# Exercice 1
mnist <- dataset_mnist()
k=5
img <- t(apply(mnist$train$x[5,,], 2, rev))
image(img,col=gray(0:255/255))
mnist$train$y[k]


# Download mnist data
mnist <- dataset_mnist()
x_train <- mnist$train$x
y_train <- mnist$train$y
x_test <- mnist$test$x
y_test <- mnist$test$y
# reshape : linearize the images
x_train <- array_reshape(x_train, c(nrow(x_train), 784))
x_test <- array_reshape(x_test, c(nrow(x_test), 784))
# rescale : rescale values from 0 to 1
x_train <- x_train / 255
x_test <- x_test / 255
# output generation
y_train <- to_categorical(y_train, 10)
y_test <- to_categorical(y_test, 10)


model=keras_model_sequential() 
#input_shape obligatoire comme premiere couche, dim. output=256
model = layer_dense(model,units = 256, activation = "relu", input_shape = 784)
#input_shape 
model = layer_dropout(model,rate = 0.4)
model = layer_dense(model,units = 128, activation = "relu") 
model = layer_dropout(model,rate = 0.3)
model = layer_dense(model,units = 10, activation = "softmax")

model=compile(model, loss = "categorical_crossentropy", optimizer = optimizer_rmsprop(), metrics = c("accuracy"))

history=fit(model, x_train, y_train, epochs = 30, batch_size = 128, validation_split = 0.2)

plot(history$metrics$acc)
points(history$metrics$val_acc,col="blue")

y_pred=predict_classes(model, x_test)

evaluate(model, x_test, y_test)
Table=table(y_pred, mnist$test$y)
sum(diag(Table))/10000



#########################################################
#########################################################
# Exercice 2

# question 1
fashion = dataset_fashion_mnist()

# question 2
f.x_train = fashion$train$x
f.y_train = fashion$train$y
f.x_test = fashion$test$x
f.y_test = fashion$test$y

rotate = function(x) t(apply(x, 2, rev))
par(mfrow=c(3,3),mgp=c(1.5,0.5,0),mar=c(3,3,2,2),cex.lab=0.9,cex.axis=0.8,cex.main=1)
for (i in 1:9) {
  image(1:28, 1:28, rotate(f.x_train[i,,]), col=grey.colors(12))
}

# question 3
# reshape
f.x_train <- array_reshape(f.x_train, c(nrow(f.x_train), 784))
f.x_test <- array_reshape(f.x_test, c(nrow(f.x_test), 784))
# rescale
f.x_train <- fa.x_train / 255
f.x_test <- fa.x_test / 255
# output generation
f.y_train <- to_categorical(fa.y_train, 10)
f.y_test <- to_categorical(fa.y_test, 10)



model=keras_model_sequential()
model = layer_dense(model,units = 200, activation = "relu", input_shape = 784)
model = layer_dropout(model,rate = 0.4)
model = layer_dense(model,units = 100, activation = "relu") 
model = layer_dropout(model,rate = 0.3)
model = layer_dense(model,units = 10, activation = "softmax")

model=compile(model, loss = "categorical_crossentropy", optimizer = optimizer_rmsprop(), metrics = c("accuracy"))

history=fit(model, x_train, y_train, epochs = 15, batch_size = 50, validation_split = 0.2)

plot(history$metrics$acc)
points(history$metrics$val_acc,col="blue")
# Accuracy final : 0.976

f.y_pred=predict_classes(model, f.x_test)
evaluate(model, f.x_test, f.y_test)
Table=table(y_pred, mnist$test$y)
sum(diag(Table))/10000
# Accuracy test : 0.976


































