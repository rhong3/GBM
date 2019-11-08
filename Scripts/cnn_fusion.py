"""
Convolutional neural nets driving code for TF2.0 for fusion

Created on 10/31/2019

@author: RH
"""
from datetime import datetime
import os
import sys
import time
import numpy as np
import tensorflow as tf
import Accessory as ac


# Define an Inception
class INCEPTION:
    # hyper parameters
    DEFAULTS = {
        "batch_size": 64,
        "dropout": 0.3,
        "learning_rate": 1E-3,
        "classes": 2
    }

    RESTORE_KEY = "cnn_to_restore"

    def __init__(self, input_dim, d_hyperparams={},
                 save_graph_def=True, meta_graph=None,
                 log_dir="./log", meta_dir="./meta", model='I1', weights = tf.constant([1., 1., 1., 1.])):

        self.input_dim = input_dim
        self.__dict__.update(INCEPTION.DEFAULTS, **d_hyperparams)
        self.sesh = tf.Session()
        self.model = model
        self.weights = weights

        if meta_graph:  # load saved graph
            model_name = os.path.basename(meta_graph)
            meta_graph = os.path.abspath(meta_graph)
            tf.train.import_meta_graph(meta_dir + '/' + model_name +'.meta').restore(
                self.sesh, meta_dir + '/' + model_name)
            handles = self.sesh.graph.get_collection(INCEPTION.RESTORE_KEY)


        else:  # build graph from scratch
            self.datetime = datetime.now().strftime(r"%y%m%d_%H%M")
            handles = self._buildGraph(self.model)
            for handle in handles:
                tf.add_to_collection(INCEPTION.RESTORE_KEY, handle)
            self.sesh.run(tf.global_variables_initializer())

        # unpack handles for tensor ops to feed or fetch for lower layers
        (self.x_in, self.is_train, self.y_in, self.logits, self.dm_in,
         self.net, self.w, self.pred, self.pred_cost, self.global_step, self.train_op, self.merged_summary) = handles

        if save_graph_def:  # tensorboard
            try:
                os.mkdir(log_dir + '/training')
                os.mkdir(log_dir + '/validation')

            except FileExistsError:
                pass

            self.train_logger = tf.summary.FileWriter(log_dir + '/training', self.sesh.graph)
            self.valid_logger = tf.summary.FileWriter(log_dir + '/validation', self.sesh.graph)

    @property
    def step(self):
        return self.global_step.eval(session=self.sesh)

    # build graph; choose a structure defined in model
    def _buildGraph(self, model):
        # image input
        x_in = tf.placeholder(tf.float32, name="x")
        x_in_reshape = tf.reshape(x_in, [-1, self.input_dim[1], self.input_dim[2], 3])
        # dropout
        dropout = self.dropout
        # label input
        y_in = tf.placeholder(dtype=tf.int32, name="y")

        # other features input
        dm_in = tf.placeholder(dtype=tf.float32, name="demographic")
        dm_in_reshape = tf.reshape(dm_in, [-1, 5])

        # train or test
        is_train = tf.placeholder_with_default(True, shape=[], name="is_train")
        classes = self.classes

        if model == 'F1':
            import IncepFusionV1
            logits, nett, ww = IncepFusionV1.incepfusionv1(x_in_reshape,
                                                           demographics=dm_in_reshape,
                                                           num_classes=classes,
                                                           is_training=is_train,
                                                           dropout_keep_prob=dropout,
                                                           scope='IncepFusionV1')
            print('Using IncepFusionV1')
        elif model == 'FS1':
            import SimpleFusionV1
            logits, nett, ww = SimpleFusionV1.simplefusionv1(x_in_reshape,
                                                           demographics=dm_in_reshape,
                                                           num_classes=classes,
                                                           is_training=is_train,
                                                           dropout_keep_prob=dropout,
                                                           scope='SimpleFusionV1')
            print('Using SimpleFusionV1')
        else:
            import IncepFusionV1
            logits, nett, ww = IncepFusionV1.incepfusionv1(x_in_reshape,
                                                           demographics=dm_in_reshape,
                                                           num_classes=classes,
                                                           is_training=is_train,
                                                           dropout_keep_prob=dropout,
                                                           scope='IncepFusionV1')
            print('Using Default: IncepFusionV1')

        pred = tf.nn.softmax(logits, name="prediction")

        global_step = tf.Variable(0, trainable=False)

        sample_weights = tf.gather(self.weights, tf.argmax(y_in, axis=1))

        pred_cost = tf.losses.softmax_cross_entropy(
            onehot_labels=y_in, logits=logits, weights=sample_weights)

        tf.summary.scalar("{}_cost".format(model), pred_cost)

        tf.summary.tensor_summary("{}_pred".format(model), pred)

        # opt = tf.optimizers.Adam(learning_rate=self.learning_rate) ---------TF2.0
        opt = tf.train.AdamOptimizer(learning_rate=self.learning_rate)
        train_op = opt.minimize(loss=pred_cost, global_step=global_step)

        merged_summary = tf.summary.merge_all()

        return (x_in, is_train,
                y_in, logits, dm_in, nett, ww, pred, pred_cost,
                global_step, train_op, merged_summary)

    # inference using trained models
    def inference(self, X, dirr, testset=None, pmd=None, train_status=False, Not_Realtest=True, bs=None):
        now = datetime.now().isoformat()[11:]
        print("------- Testing begin: {} -------\n".format(now))
        rd = 0
        pdx = []
        yl = []
        if Not_Realtest:
            itr, file, ph = X.data(train=False)
            next_element = itr.get_next()
            with tf.Session() as sessa:
                sessa.run(itr.initializer, feed_dict={ph: file})
                while True:
                    try:
                        x, y, dm = sessa.run(next_element)
                        feed_dict = {self.x_in: x, self.demographic: dm, self.is_train: train_status}
                        fetches = [self.pred, self.net, self.w]
                        pred, net, w = self.sesh.run(fetches, feed_dict)
                        # ac.CAM(net, w, pred, x, y, dirr, 'Test', bs, pmd, rd)
                        net = np.mean(net, axis=(1, 2))
                        if rd == 0:
                            pdx = pred
                            yl = y
                            netl = net
                        else:
                            pdx = np.concatenate((pdx, pred), axis=0)
                            yl = np.concatenate((yl, y), axis=0)
                            netl = np.concatenate((netl, net), axis=0)
                        rd += 1
                    except tf.errors.OutOfRangeError:
                        ac.metrics(pdx, yl, dirr, 'Test', pmd, testset)
                        ac.tSNE_prep(flatnet=netl, ori_test=testset, y=yl, pred=pdx, path=dirr, pmd=pmd)
                        break
        else:
            itr, file, ph = X.data(Not_Realtest=False, train=False)
            next_element = itr.get_next()
            with tf.Session() as sessa:
                sessa.run(itr.initializer, feed_dict={ph: file})
                while True:
                    try:
                        x, dm  = sessa.run(next_element)
                        feed_dict = {self.x_in: x, self.demographic: dm, self.is_train: train_status}
                        fetches = [self.pred, self.net, self.w]
                        pred, net, w = self.sesh.run(fetches, feed_dict)
                        ac.CAM_R(net, w, pred, x, dirr, 'Test', bs, rd)
                        if rd == 0:
                            pdx = pred
                        else:
                            pdx = np.concatenate((pdx, pred), axis=0)
                        rd += 1
                    except tf.errors.OutOfRangeError:
                        ac.realout(pdx, dirr, 'Test', pmd)
                        break

        now = datetime.now().isoformat()[11:]
        print("------- Testing end: {} -------\n".format(now), flush=True)

    # training
    def train(self, X, VAX, ct, bs, dirr, pmd, max_iter=np.inf,
              cross_validate=True, save=True, outdir="./out"):
        start_time = time.time()
        svs = 0
        if save:
            saver = tf.train.Saver(tf.global_variables(), max_to_keep=None)

        try:
            err_train = 0
            now = datetime.now().isoformat()[11:]
            print("------- Training begin: {} -------\n".format(now))
            itr, file, ph = X.data()
            next_element = itr.get_next()

            vaitr, vafile, vaph = VAX.data(train=False)
            vanext_element = vaitr.get_next()

            with tf.Session() as sessa:
                sessa.run(itr.initializer, feed_dict={ph: file})
                sessa.run(vaitr.initializer, feed_dict={vaph: vafile})
                train_cost = []
                validation_cost = []
                valid_cost = 0
                while True:
                    try:
                        x, y, dm = sessa.run(next_element)
                        feed_dict = {self.x_in: x, self.demographic: dm}

                        fetches = [self.merged_summary, self.logits, self.pred,
                                   self.pred_cost, self.global_step, self.train_op]

                        summary, logits, pred, cost, i, _ = self.sesh.run(fetches, feed_dict)

                        self.train_logger.add_summary(summary, i)
                        err_train += cost

                        if i < 2:
                            train_cost.append(cost)

                        try:
                            mintrain = min(train_cost)
                        except ValueError:
                            mintrain = 0

                        if cost <= mintrain and i > 29999:
                            temp_valid = []
                            for iii in range(100):
                                x, y, dm = sessa.run(vanext_element)
                                feed_dict = {self.x_in: x, self.demographic: dm, self.is_train: False}
                                fetches = [self.pred_cost, self.merged_summary]
                                valid_cost, valid_summary = self.sesh.run(fetches, feed_dict)
                                self.valid_logger.add_summary(valid_summary, i)
                                temp_valid.append(valid_cost)

                            tempminvalid = np.mean(temp_valid)
                            try:
                                minvalid = min(validation_cost)
                            except ValueError:
                                minvalid = 0

                            if tempminvalid <= minvalid:
                                train_cost.append(cost)
                                print("round {} --> loss: ".format(i), cost)
                                print("round {} --> validation loss: ".format(i), tempminvalid)
                                print("New Min loss model found!", flush=True)
                                validation_cost.append(tempminvalid)
                                if save:
                                    outfile = os.path.join(os.path.abspath(outdir),
                                                           "{}_{}".format(self.model,
                                                                          "_".join(['dropout', str(self.dropout)])))
                                    saver.save(self.sesh, outfile, global_step=None)
                                    svs = i

                        else:
                            train_cost.append(cost)

                        if i % 1000 == 0:
                            print("round {} --> loss: ".format(i), cost)
                            temp_valid = []
                            for iii in range(100):
                                x, y, dm = sessa.run(vanext_element)
                                feed_dict = {self.x_in: x, self.demographic: dm, self.is_train: False}
                                fetches = [self.pred_cost, self.merged_summary]
                                valid_cost, valid_summary = self.sesh.run(fetches, feed_dict)
                                self.valid_logger.add_summary(valid_summary, i)
                                temp_valid.append(valid_cost)
                            tempminvalid = np.mean(temp_valid)
                            try:
                                minvalid = min(validation_cost)
                            except ValueError:
                                minvalid = 0
                            validation_cost.append(tempminvalid)
                            print("round {} --> Step Average validation loss: ".format(i), tempminvalid, flush=True)

                            if save and tempminvalid <= minvalid:
                                print("New Min loss model found!")
                                print("round {} --> loss: ".format(i), cost)
                                outfile = os.path.join(os.path.abspath(outdir),
                                                       "{}_{}".format(self.model,
                                                                      "_".join(['dropout', str(self.dropout)])))
                                saver.save(self.sesh, outfile, global_step=None)
                                svs = i

                            if i > 99999:
                                valid_mean_cost = np.mean(validation_cost[-10:-1])
                                print('Mean validation loss: {}'.format(valid_mean_cost))
                                if valid_cost > valid_mean_cost:
                                    print("Early stopped! No improvement for at least 10000 iterations", flush=True)
                                    break
                                else:
                                    print("Passed early stopping evaluation. Continue training!")

                        if i >= max_iter-2:
                            print("final avg loss (@ step {} = epoch {}): {}".format(
                                i + 1, np.around(i / ct * bs), err_train / i))

                            now = datetime.now().isoformat()[11:]
                            print("------- Training end: {} -------\n".format(now))

                            now = datetime.now().isoformat()[11:]
                            print("------- Final Validation begin: {} -------\n".format(now))
                            x, y, dm = sessa.run(vanext_element)
                            feed_dict = {self.x_in: x, self.demographic: dm, self.is_train: False}
                            fetches = [self.pred_cost, self.merged_summary]
                            valid_cost, valid_summary= self.sesh.run(fetches, feed_dict)

                            self.valid_logger.add_summary(valid_summary, i)
                            print("round {} --> Final Last validation loss: ".format(i), valid_cost)
                            now = datetime.now().isoformat()[11:]
                            print("------- Final Validation end: {} -------\n".format(now), flush=True)
                            try:
                                self.train_logger.flush()
                                self.train_logger.close()
                                self.valid_logger.flush()
                                self.valid_logger.close()

                            except AttributeError:  # not logging
                                print('Not logging')
                            break

                    except tf.errors.OutOfRangeError:
                        print("final avg loss (@ step {} = epoch {}): {}".format(
                            i + 1, np.around(i / ct * bs), err_train / i))

                        now = datetime.now().isoformat()[11:]
                        print("------- Training end: {} -------\n".format(now))

                        now = datetime.now().isoformat()[11:]
                        print("------- Final Validation begin: {} -------\n".format(now))
                        x, y, dm = sessa.run(vanext_element)
                        feed_dict = {self.x_in: x, self.demographic: dm, self.is_train: False}
                        fetches = [self.pred_cost, self.merged_summary, self.pred, self.net, self.w]
                        valid_cost, valid_summary, pred, net, w = self.sesh.run(fetches, feed_dict)

                        self.valid_logger.add_summary(valid_summary, i)
                        print("round {} --> Final Last validation loss: ".format(i), valid_cost)
                        ac.CAM(net, w, pred, x, y, dirr, 'Validation', bs, pmd)
                        ac.metrics(pred, y, dirr, 'Validation', pmd)
                        now = datetime.now().isoformat()[11:]
                        print("------- Final Validation end: {} -------\n".format(now), flush=True)

                        try:
                            self.train_logger.flush()
                            self.train_logger.close()
                            self.valid_logger.flush()
                            self.valid_logger.close()

                        except AttributeError:  # not logging
                            print('Not logging')

                        break
                try:
                    print("final avg loss (@ step {} = epoch {}): {}".format(
                        i + 1, np.around(i / ct * bs), err_train / i))

                    now = datetime.now().isoformat()[11:]
                    print("------- Training end: {} -------\n".format(now))

                    if svs < 15000 and save:
                            print("Save the last model as the best model.")
                            outfile = os.path.join(os.path.abspath(outdir),
                                                   "{}_{}".format(self.model, "_".join(['dropout', str(self.dropout)])))
                            saver.save(self.sesh, outfile, global_step=None)

                    now = datetime.now().isoformat()[11:]
                    print("------- Validation begin: {} -------\n".format(now))
                    x, y, dm = sessa.run(vanext_element)
                    feed_dict = {self.x_in: x, self.demographic: dm, self.is_train: False}
                    fetches = [self.pred_cost, self.merged_summary, self.pred, self.net, self.w]
                    valid_cost, valid_summary, pred, net, w = self.sesh.run(fetches, feed_dict)

                    self.valid_logger.add_summary(valid_summary, i)
                    print("round {} --> Last validation loss: ".format(i), valid_cost)
                    ac.CAM(net, w, pred, x, y, dirr, 'Validation', bs, pmd)
                    ac.metrics(pred, y, dirr, 'Validation', pmd)
                    now = datetime.now().isoformat()[11:]
                    print("------- Validation end: {} -------\n".format(now), flush=True)

                    try:
                        self.train_logger.flush()
                        self.train_logger.close()
                        self.valid_logger.flush()
                        self.valid_logger.close()

                    except AttributeError:  # not logging
                        print('Not logging')

                except tf.errors.OutOfRangeError:
                    print("final avg loss (@ step {} = epoch {}): {}".format(
                        i + 1, np.around(i / ct * bs), err_train / i))

                    now = datetime.now().isoformat()[11:]
                    print("------- Training end: {} -------\n".format(now))
                    print('No more validation needed!')

            print("--- %s seconds ---" % (time.time() - start_time))

        except KeyboardInterrupt:

            print("final avg loss (@ step {} = epoch {}): {}".format(
                i, np.around(i / ct * bs), err_train / i))

            now = datetime.now().isoformat()[11:]
            print("------- Training end: {} -------\n".format(now))

            if save:
                outfile = os.path.join(os.path.abspath(outdir),
                                       "{}_{}".format(self.model, "_".join(['dropout', str(self.dropout)])))
                saver.save(self.sesh, outfile, global_step=None)
            try:
                self.train_logger.flush()
                self.train_logger.close()
                self.valid_logger.flush()
                self.valid_logger.close()

            except AttributeError:  # not logging
                print('Not logging')

            sys.exit(0)











