import matplotlib
# matplotlib.use('Agg')  # no UI backend

import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt


mnist = tf.keras.datasets.mnist  # 28x28 px images of handwritten digits 0-9

(x_train, y_train), (x_test, y_test) = mnist.load_data()
x_train = tf.keras.utils.normalize(x_train, axis=1)
x_test = tf.keras.utils.normalize(x_test, axis=1)

# displays the first image in the data set
# plt.imshow(x_train[0], cmap=plt.cm.binary)
# plt.show()

model = tf.keras.models.Sequential()  # feed-forward model
model.add(tf.keras.layers.Flatten())  # input layer, flattens data
model.add(tf.keras.layers.Dense(128, activation=tf.nn.relu))  # 128-neuron layer, with rectified-linear activation function
model.add(tf.keras.layers.Dense(128, activation=tf.nn.relu))  # third layer, same structure as the second
model.add(tf.keras.layers.Dense(10, activation=tf.nn.softmax))  # output layer

model.compile(
    optimizer='adam',  # method for minimizing the los
    loss='sparse_categorical_crossentropy',  # method for calculating loss
    metrics=['accuracy']
)

model.fit(x_train, y_train, epochs=3)

# validation
val_loss, val_acc = model.evaluate(x_test, y_test)
print(val_loss)
print(val_acc)

# save a model
model.save('number_reader.model')

# load a model
new_model = tf.keras.models.load_model('number_reader.model')

# make a prediction
predictions = new_model.predict(x_test)
print(np.argmax(predictions[0]))
plt.imshow(x_test[0], cmap=plt.cm.binary)
plt.show()
