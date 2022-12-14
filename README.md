# Extraction of collagen morphological features
The program is used to extract the morphological features of collagen in multiphoton images, including:
* Collagen proportionate area
* Collagen fiber number
* Collagen fiber length
* Collagen fiber width
* Collagen fiber straightness
* Collagen fiber cross-link density
* Collagen fiber cross-link space
* Collagen fiber orientation

## Dependency and installation
* MATLAB 2016b   

Firstly, download the whole folder, in which the photo folder is the place where the image is stored, and the Mask folder is used to store the output of collagen after the program is segmented. The main program is ```main```, where two parameters are input according to the actual images.

Image brightness threshold:  
```
threshold = **; # The range is 0-255, and the appropriate threshold is selected based on the brightness of the image     
```
Resolution ratio:   
```
MMP = ** ; # Resolution ratio, which represents the actual size of a pixel in micrometers. 
```
After setting the parameters, ```run``` matlab software, and the extracted values will be stored in the ```FEA``` of the parameter list in order

## Citation

```
@inproceedings{Feature extraction program,
  title={Rapid and label-free detection of gastrointestinal stromal tumor via a combination of two-photon microscopy and imaging analysis},
  author={Lianhuang Li, Xingxin Huang, Shichao Zhang, Zhenlin Zhan, Deyong Kang, Guoxian Guan, Shuoyu Xu, Yongjian Zhou, Jianxin Chen },
  booktitle={BMC},
  year={2022}
}
```

