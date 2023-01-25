import numpy as np
import tensorflow as tf
from tensorflow import keras
import sys
import datetime
import os

if "snakemake" not in globals():
    class Snakemake:
        input = [sys.argv[1]]
        output = [sys.argv[2]]
        params = [sys.argv[2]]
    snakemake = Snakemake()

data = tf.data.Dataset.load(os.path.dirname(snakemake.input[0]))#"tokenized_5mer.tf")#data_path)
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))
tf.debugging.set_log_device_placement(True)
vocab_sz = 5590725 #100000#14162 # 4194304 possible kmers, vocab size = # 11mers seen in vocab file
embedding_dim = 45 
max_len = 31 
bsz = 10**4 # batch size
model = keras.Sequential()
model.add(keras.layers.Embedding(vocab_sz, embedding_dim, input_length=max_len, trainable=True))
#model.add(keras.layers.Conv1D(17, 3, activation="relu"))
model.add(keras.layers.Bidirectional(keras.layers.LSTM(128)))
model.add(keras.layers.Dense(48,activation="relu"))
model.add(keras.layers.Dropout(0.1))
model.add(keras.layers.Dense(16, activation="tanh"))
model.add(keras.layers.Dense(1, activation="tanh"))
model.compile(loss="mean_squared_error",optimizer="adam",metrics=["accuracy"])

#data_b = tf.data.Dataset.load("tokenized_refined_bad.tf")
sz = len(data) #100#len(data_g) + len(data_b)
print(sz)
train_sz = int(sz*0.8) #int(len(data)*0.8) # 80-10-10 train-validation-test split
print(train_sz)
test_sz = int(sz*0.1)

#data = tf.data.Dataset.sample_from_datasets([data_g, data_b], weights=[0.5, 0.5]) #data.take(train_sz).batch(1000)
#print(data)
data = data.shuffle(10000, seed=42)
train = data.take(train_sz).batch(bsz)
val = data.skip(train_sz).take(test_sz).batch(bsz)
#val.save(datetime.datetime.now().strftime("m-d-H-M")+".val.tf")
test = data.skip(train_sz).skip(test_sz).take(test_sz).batch(bsz)
del data

nepochs = 15
model.fit(train, epochs=nepochs, validation_data=val)

model_name = snakemake.params[0]
model.save(model_name)

with open(model_name + ".val.loss","w") as f:
    f.write(str(model.history.history))

with open(model_name + ".unseen.loss","w") as f:
    f.write(str(model.evaluate(test)))

seen_test = train.unbatch().shuffle(1000,seed=42).take(test_sz).batch(bsz)
with open(model_name +".seen.loss","w") as f:
    f.write(str(model.evaluate(seen_test)))

