import numpy as np
from scipy.io import loadmat
from matplotlib import pyplot as plt


def showmat(path):
    img = loadmat(path)

    # TODO: there has to be a better way of getting the array out of the loadmat return object
    keys = img.keys()
    key = None
    for k in keys:
        if k[0] != '_':
            key = k
            break
    if key is None:
        raise Exception("Error loading .mat file.")

    imgplot = plt.imshow(img[key][:100,:100], cmap='gray')
    plt.show()


def shownpy(path):
    img = np.load(path)
    imgplot = plt.imshow(img, cmap='gray')
    plt.show()


def test():
    pathmat = r'./DeepMCR_package_v0.1/imgs/sandstone.mat'
    showmat(pathmat)
    for i in range(1, 136, 10):
        pathnpy = './DeepMCR_package_v0.1/.output/record_step_' + str(i) + '00.npy'
        shownpy(pathnpy)
    showmat(pathmat)

if __name__ == '__main__':
    # TODO: organize testing functionality a bit better
    test()
