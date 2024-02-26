import tensorflow.compat.v1 as tf1
import tensorflow as tf



class DeepCas9(object):

    def __init__(self, filter_size, filter_num, node_1=80, node_2=60, l_rate=0.005):
        # filter_size = [3, 5, 7]
        # filter_num = [100, 70, 40]
        length = 30
        self.inputs = tf1.placeholder(tf1.float32, [None, 1, length, 4])
        self.targets = tf1.placeholder(tf1.float32, [None, 1])
        self.is_training = tf1.placeholder(tf1.bool)

        def create_new_conv_layer(input_data, num_input_channels, num_filters, filter_shape, pool_shape, name):
            # setup the filter input shape for tf.nn.conv_2d
            conv_filt_shape = [filter_shape[0], filter_shape[1], num_input_channels,
                               num_filters]

            # initialise weights and bias for the filter
            weights = tf1.Variable(tf1.truncated_normal(conv_filt_shape, stddev=0.03),
                                   name=name + '_W')
            bias = tf1.Variable(tf1.truncated_normal([num_filters]), name=name + '_b')

            # setup the convolutional layer operation
            out_layer = tf.nn.conv2d(input_data, weights, [1, 1, 1, 1], padding='VALID')
            #out_layer = tf1.nn.conv2d(input_data, weights, [1, 1, 1, 1], padding='VALID')

            # add the bias
            out_layer += bias

            # apply a ReLU non-linear activation
            #out_layer = tf1.layers.dropout(tf1.nn.relu(out_layer), 0.3, self.is_training)
            out_layer = tf.keras.layers.Dropout(rate=0.3)(tf.nn.relu(out_layer))

            # now perform max pooling
            ksize = [1, pool_shape[0], pool_shape[1], 1]
            strides = [1, 1, 2, 1]
            out_layer = tf.nn.avg_pool(out_layer, ksize=ksize, strides=strides,
                                        padding='SAME')
            return out_layer

        # def end: create_new_conv_layer

        L_pool_0 = create_new_conv_layer(self.inputs, 4, filter_num[0], [1, filter_size[0]], [1, 2], name='conv1')
        L_pool_1 = create_new_conv_layer(self.inputs, 4, filter_num[1], [1, filter_size[1]], [1, 2], name='conv2')
        L_pool_2 = create_new_conv_layer(self.inputs, 4, filter_num[2], [1, filter_size[2]], [1, 2], name='conv3')

        with tf1.variable_scope('Fully_Connected_Layer1'):
            layer_node_0 = int((length - filter_size[0]) / 2) + 1
            node_num_0 = layer_node_0 * filter_num[0]
            layer_node_1 = int((length - filter_size[1]) / 2) + 1
            node_num_1 = layer_node_1 * filter_num[1]
            layer_node_2 = int((length - filter_size[2]) / 2) + 1
            node_num_2 = layer_node_2 * filter_num[2]

            L_flatten_0 = tf.reshape(L_pool_0, [-1, node_num_0])
            L_flatten_1 = tf.reshape(L_pool_1, [-1, node_num_1])
            L_flatten_2 = tf.reshape(L_pool_2, [-1, node_num_2])
            L_flatten = tf.concat([L_flatten_0, L_flatten_1, L_flatten_2], 1, name='concat')

            node_num = node_num_0 + node_num_1 + node_num_2
            W_fcl1 = tf1.get_variable("W_fcl1", shape=[node_num, node_1])
            B_fcl1 = tf1.get_variable("B_fcl1", shape=[node_1])
            L_fcl1_pre   = tf.nn.bias_add(tf.matmul(L_flatten, W_fcl1), B_fcl1)
            L_fcl1       = tf.nn.relu(L_fcl1_pre)
            L_fcl1_drop  = tf.keras.layers.Dropout(rate=0.3)(L_fcl1)


        with tf1.variable_scope('Fully_Connected_Layer2'):
            W_fcl2 = tf1.get_variable("W_fcl2", shape=[node_1, node_2])
            B_fcl2 = tf1.get_variable("B_fcl2", shape=[node_2])
            L_fcl2_pre = tf.nn.bias_add(tf.matmul(L_fcl1_drop, W_fcl2), B_fcl2)
            L_fcl2 = tf.nn.relu(L_fcl2_pre)
            L_fcl2_drop = tf.keras.layers.Dropout(rate=0.3)(L_fcl2)


        with tf1.variable_scope('Output_Layer'):
            W_out = tf1.get_variable("W_out",
                                     shape=[node_2, 1])  # , initializer=tf.contrib.layers.xavier_initializer())
            B_out = tf1.get_variable("B_out", shape=[1])  # , initializer=tf.contrib.layers.xavier_initializer())
            self.outputs = tf.nn.bias_add(tf.matmul(L_fcl2_drop, W_out), B_out)

        # Define loss function and optimizer
        self.obj_loss = tf.reduce_mean(tf.square(self.targets - self.outputs))
        self.optimizer = tf1.train.AdamOptimizer(l_rate).minimize(self.obj_loss)
    # def end: def __init__
