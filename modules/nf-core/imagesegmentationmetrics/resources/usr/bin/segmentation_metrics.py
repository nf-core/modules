#!/usr/bin/env python3

import monai.metrics as mmetrics
import tifffile
import torch
import numpy as np
import pandas as pd
from enum import Enum
import fire
from time import time
import seg_metrics.seg_metrics as sg
import logging

VERSION = "0.1.1"

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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
    return torch.reshape(torch.from_numpy(im), [1] + [1] + list(im.shape)).float()


def load_dataset(path_ref, path_list_seg) -> tuple[torch.Tensor, torch.Tensor, np.ndarray, np.ndarray]:
    """
    Load the dataset and return a list of torch tensors.
    """
    imref = read_image(path_ref)
    imarray = []
    for path_seg in path_list_seg.split(' '):
        imarray.append(read_image(path_seg))
    imref_tor = torch.cat([convert_to_tensor(binarize(imref))] * len(imarray), 0)
    imarray_tor = torch.cat([convert_to_tensor(binarize(im)) for im in imarray], 0)
    imref_np = binarize(imref)
    imarray_np = [binarize(im) for im in imarray]
    return imref_tor, imarray_tor, imref_np, imarray_np

def parse_metrics(str_metrics) -> list[tuple[str, str]]:
    """
    Parse the metrics from the list of metrics.
    """
    if str_metrics == "":
        return [('monai', i) for i in monai_metrics.__members__] + [('sgm', i) for i in sgm_metrics]
    list_metrics = []
    for i in str_metrics.split(' '):
        if i == "monai":
            list_metrics += [(i, j) for j in monai_metrics.__members__]
        elif i == "sgm":    
            list_metrics += [(i, j) for j in sgm_metrics]
        else:
            pkg_met = i.split('_')
            if len(pkg_met) == 1:
                list_metrics.append(('monai',pkg_met[0])) # Default to monai
            else:
                list_metrics.append((pkg_met[0], pkg_met[1]))
    return list_metrics

def benchmark_metrics(ref_tor, seg_tor, ref_np, seg_np, list_metrics) -> list[torch.Tensor]:
    """
    Compute the benchmark metrics.
    Add RAM/CPU benchmark
    """
    list_results = []
    for i,j in list_metrics:
        ts = time() # Start timing the operation
        if i == 'monai':
            list_results.append(getattr(monai_metrics, j)(ref_tor,seg_tor))
        elif i == 'sgm':
            list_results.append(get_sgm_metrics(ref_np, seg_np, j))
        logger.info(f"{j}({i}) runtime: {(time()-ts)} sec")
    return(list_results)

def get_results(results, path_list_seg, list_metrics) -> pd.DataFrame:
    """
    Get the results of the metrics.
    """
    out = pd.DataFrame(torch.cat(results, 1).numpy(),
        index = path_list_seg.split(' '), columns = map(lambda x: f"{x[0]}_{x[1]}", list_metrics))
    out.to_csv('segmentation_benchmark_results.csv')
    return out

def get_sgm_metrics(ref_np, seg_np_list, metric_name) -> torch.Tensor:
    """
    Getting score for seg-metrics library metrics
    """
    tensor_res = torch.empty(len(seg_np_list),1)
    for i,seg_np in enumerate(seg_np_list):
        metrics = sg.write_metrics(labels=[1],  # exclude background if needed
                      gdth_img=ref_np,
                      pred_img=seg_np,
                      metrics=[metric_name])
        tensor_res[i][0] = metrics[0][metric_name][0]
    return tensor_res
    
sgm_metrics = ['dice', 'jaccard', 'precision', 'recall', 'fpr', 'fnr', 'hd', 'msd', 'stdsd']

class monai_metrics(Enum):
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
    MSEMetric = mmetrics.MSEMetric()
    MAEMetric = mmetrics.MAEMetric()
    RMSEMetric = mmetrics.RMSEMetric()
    PSNRMetric = mmetrics.PSNRMetric(max_val = 1)
    

def main(path_ref, path_list_seg, str_metrics) -> None:
    """
    Main function to compute the metrics.
    input:
    path_ref: path to the ground truth segmentation (as TIFF file)
    path_list_seg: space separated list of paths to the segmentations to compare
    list_metrics: space separated list of metrics to compute
    current metrics: monai_DiceMetric monai_MeanIoU monai_GeneralizedDiceScore monai_HausdorffDistanceMetric monai_SurfaceDistanceMetric monai_SurfaceDiceMetric monai_MSEMetric monai_RMSEMetric monai_PSNRMetric sgm_dice sgm_jaccard sgm_precision sgm_recall sgm_fpr sgm_fnr sgm_hd sgm_msd sgm_stdsd
    """
    ref_tor, seg_tor, ref_np, seg_np = load_dataset(path_ref, path_list_seg)
    list_metrics = parse_metrics(str_metrics)
    results = benchmark_metrics(ref_tor, seg_tor, ref_np, seg_np, list_metrics)
    get_results(results, path_list_seg, list_metrics)

if __name__ == '__main__':
    options = {
        "run" : main,
        "version" : VERSION
    }
    fire.Fire(options)