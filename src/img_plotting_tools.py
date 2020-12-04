import numpy as np
from scipy.io import loadmat
from matplotlib import pyplot as plt


def showmat(path):
    img = loadmat(path)
    keys = img.keys()
    key = None
    for k in keys:
        if k[0] != '_':
            key = k
            break
    if key is None:
        raise Exception("Error loading .mat file.")
    imgplot = plt.imshow(img[key], cmap='gray')
    plt.show()


def shownpy(path):
    img = np.load(path)
    imgplot = plt.imshow(img, cmap='gray')
    plt.show()
