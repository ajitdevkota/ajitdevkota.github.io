---
title: Image Classification using a Convolutional Neural Network
date: 2021-05-15 14:29:20 +/-0000
categories: [Blogging, Demonstration]
tags: [tensorflow, jupyter, cnn, cifar-10]
---

## Overview
This project aims to design and train a convolutional neural network (CNN) on the CIFAR-10
dataset. Model over-fitting and poor performance are the main issues in designing a CNN.
To solve these problems, I conducted a list of well-thought strategies.

![Window shadow](/assets/img/posts/cnn-cifar-10.png){: .shadow style="max-width: 90%" }
_CIFAR-10 Dataset_

**Workflows used**: Git and Agile.
**Technologies used**: Python, Jupyter, TensorFlow.

**Note**: *This project was done as part of assessment for module CSC2034.*

## Task
Design and train a Convolutional Neural Network that performs as well as possible on the CIFAR-10 dataset.

## Results
After implementing pre-defined CNN architecture, optimisation algorithm, image data augmentation techniques and tuning general model parameters, overall model train accuracy jumped to 94%, and test accuracy reached 89%. Both train and test losses were below 0.4.

![Window shadow](/assets/img/posts/cnn-final-model-accuracy-loss.jpg){: .shadow style="max-width: 90%" }
_Accuracy and loss of the final model_

![Window shadow](/assets/img/posts/cnn-final-model-confusion-matrix.jpg){: .shadow style="max-width: 90%" }
_Confusion matrix of the final model_

## Summary of CNN Development
In this study, I performed systematic experiments to design and train a convolutional neural network.
Throughout the development of the model, I improved CNN performance and reduced overfitting with
four different strategies: by selecting pre-defined CNN architecture, optimisation algorithm,
implementing image data augmentation and tuning general training hyper-parameters such as batch
size and the number of epochs.

![Window shadow](/assets/img/posts/cnn-development-summary.png){: .shadow style="max-width: 90%" }
_Summary of CNN development_

## Learning Outcomes
- Selected and applied techniques appropriate to a particular task or goal.
- Demonstrated use of relevant software tools.
- Critically evaluated and assessed research papers.
- Gained knowledge of the basic concepts of deep learning.

## Learn More
For more knowledge about the project, visit my project's [repository](https://github.com/mantas2000/CIFAR-10-CNN){: target="_blank"}.
