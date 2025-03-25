#!/usr/bin/env python3

import monai.metrics as mmetrics
import tifffile
import torch
import numpy as np
import pandas as pd
from enum import Enum
import fire
from time import time

VERSION = "0.0.1"

def read_image(path) -> np.ndarray:
    """
    Read a image from a given path.
    For now only tiff files are supported.
    """
    return tifffile.imread(path)

def binarize(im) -> np.ndarray:
    """
    Image conversion
    For now only binarize the image.
    """
    return (im > 0).astype(int)

def convert_to_tensor(im) -> torch.Tensor:
    """
    Convert a numpy array to a torch tensor with padding for batch and channel.
    """
    return torch.reshape(torch.from_numpy(im), [1] + [1] + list(im.shape))


def load_dataset(path_ref, path_list_seg) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Load the dataset and return a list of torch tensors.
    """
    imref = read_image(path_ref)
    imarray = []
    for path_seg in path_list_seg.split(' '):
        imarray.append(read_image(path_seg))
    imref = torch.cat([convert_to_tensor(binarize(imref))] * len(imarray), 0)
    imarray = torch.cat([convert_to_tensor(binarize(im)) for im in imarray], 0)
    return imref, imarray

def benchmark_metrics(ref, seg, list_metrics) -> list[torch.Tensor]:
    """
    Compute the benchmark metrics.
    Add RAM/CPU benchmark
    """
    list_results = []
    for i in list_metrics.split(' '):
        ts = time() # Start timing the operation
        list_results.append(getattr(Metrics, i)(ref,seg))
        print(f"{i} runtime: {(time()-ts)} sec")
    return(list_results)

def get_results(results, path_list_seg, list_metrics) -> pd.DataFrame:
    """
    Get the results of the metrics.
    """
    out = pd.DataFrame(torch.cat(results, 1).numpy(),
        index = path_list_seg.split(' '), columns = list_metrics.split(' '))
    out.to_csv('segmentation_benchmark_results.csv')
    return out


class Metrics(Enum):
    """
    Enum class for the metrics.
    Add function to change metrics specific parameters here.
    """
    def __call__(self, *args, **kwargs):
        return self.value(*args, **kwargs)
    DiceMetric = mmetrics.DiceMetric()
    MeanIoU = mmetrics.MeanIoU()
    GeneralizedDiceScore = mmetrics.GeneralizedDiceScore(reduction='mean')
    HausdorffDistanceMetric = mmetrics.HausdorffDistanceMetric()
    SurfaceDistanceMetric = mmetrics.SurfaceDistanceMetric()
    SurfaceDiceMetric = mmetrics.SurfaceDiceMetric(class_thresholds = [2])
    # PanopticQualityMetric = mmetrics.PanopticQualityMetric(num_classes = 2)
    MSEMetric = mmetrics.MSEMetric()
    # MAEMetric = mmetrics.MAEMetric()
    RMSEMetric = mmetrics.RMSEMetric()
    PSNRMetric = mmetrics.PSNRMetric(max_val = 1)

def main(path_ref, path_list_seg, list_metrics):
    """
    Main function to compute the metrics.
    input:
    path_ref: path to the ground truth segmentation (as TIFF file)
    path_list_seg: space separated list of paths to the segmentations to compare
    list_metrics: space separated list of metrics to compute
    current metrics: DiceMetric MeanIoU GeneralizedDiceScore HausdorffDistanceMetric SurfaceDistanceMetric SurfaceDiceMetric MSEMetric RMSEMetric PSNRMetric
    """
    ref, seg = load_dataset(path_ref, path_list_seg)
    results = benchmark_metrics(ref, seg, list_metrics)
    get_results(results, path_list_seg, list_metrics)

if __name__ == '__main__':
    options = {
        "run" : main,
        "version" : VERSION
    }
    fire.Fire(options)
